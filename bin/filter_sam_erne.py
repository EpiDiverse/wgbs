#!/usr/bin/env python

'''
Title: filter_sam_erne.py
Date: 20180918
Author: Adam Nunn
Description:
	This program takes an input BAM file, typically from erne-bs5, and filters the
	alignments to exclude beyond-end-of-reference alignments and those that exceed
	a user-defined threshold of opposite-strand mismatches. Write all other
	alignments to a new BAM file. Input BAM file must be sorted and indexed for
	parallel processing.

List of functions:
	main()
	build_genome()
	worker()
	merger()
	recursive_merger()

Procedure:
	1. Build FASTA dictionary
	2. Declare variables and build threads pool
	3. Open initial BAM file instance to read references
	4. Fire off workers to process reads for each scaffold 'ref' from 'BAM'
	5. Merge results from the workers with a recursive merging function

Usage:
		./filter_sam_erne.py [-h, --help] \
			[-c, --cutoff] float \
			[-t, --threads] int \
			[-T, --temp] <TEMP> \
			[infasta] [infile] [outfile] 
	eg. ./filter_sam_erne.py -c 0.05 -t 4 -T /var/tmp in.fasta in.bam out.bam
'''

###################
## INIT ENVIRONMENT

import multiprocessing as mp
import resource, tempfile, os, argparse
import pysam


##################
## DEFINE __MAIN__
def main(FASTA,BAM,OUT,TEMP,CUTOFF,THREADS):

	# 1) Build FASTA dictionary
	genome = build_genome(FASTA)

	# 2) Declare variables and build threads pool
	pool = mp.Pool(int(THREADS))
	jobs = []

	# determine open file limit
	rlimit = resource.getrlimit(resource.RLIMIT_NOFILE)
	rlimit = rlimit[0]

	# define temp directory
	TDIR = tempfile.mkdtemp(dir=TEMP)

	# 3) Open initial BAM file instance to read references
	with pysam.AlignmentFile(BAM, "rb") as original:

		header = original.header.to_dict()
		if "CO" in header: del header["CO"]

		lengths = original.lengths
		for ref in original.get_index_statistics():

			# 4) Fire off workers to process reads for each scaffold 'ref' from 'BAM'
			if ref.mapped > 0:
				tid = original.get_tid(ref.contig)
				job = pool.apply_async(worker, (BAM, TDIR, CUTOFF, header, ref.contig, genome.pop(ref.contig), lengths[tid]))
				jobs.append(job)
	
	# 5) Merge results from the workers with a recursive merging function
	bam = recursive_merger(jobs,rlimit,pool,TDIR)
	bam = bam[0].get()
	os.replace(bam, OUT)

	# close the pool
	pool.close()
	pool.join()


## END OF __MAIN__
##################

###################
## DEFINE FUNCTIONS

####### Function to build 'genome' dictionary from 'FASTA' file object
def build_genome(FASTA):

	# Open FASTA file object for reading
	with open(FASTA, 'r') as fasta:

		genome = dict()
		first = True

		# Iterate through reference FASTA to build genome dictionary
		for line in fasta:
			line = line.rstrip()

			# Get the first sequence ID
			if line.startswith(">") and (first == True):
				line = line.split(" ")
				chr = line[0][1:]
				longline = ''
				first = False

			# Get all other sequence IDs
			elif line.startswith(">") and (first == False):
				genome[chr] = longline.upper()
				line = line.split(" ")
				chr = line[0][1:]
				longline = ''

			# Concatenate the sequence lines
			else: longline += line

		# Add the final sequence entry to the dictionary
		genome[chr] = longline.upper()
	
	return genome



####### Function for READING the input reads, modifying, and sending them to 'q'
def worker(BAM,TDIR,CUTOFF,header,rnam,rseq,rlen):

	'''
	BAM = path to input bam file eg. "/path/to/input.bam"
	TDIR = path to temp dir for processed bam files eg. "/var/tmp/adfas7d"
	CUTOFF = float proportion eg. 0.03
	rnam = current scaffold or chromosome reference eg. "Chr1"
	rseq = reference sequence result of genome[rnam] eg. "CACGTACGTAGCTAGCTAAA..."
	rlen = reference sequence length integer eg. 12512
	'''

	# generate a name for the modified bam file
	name = TDIR + "/" + rnam + ".bam"

	# Open pysam.AlignmentFile objects for reading and writing
	with pysam.AlignmentFile(BAM, "rb") as original, pysam.AlignmentFile(name, "wb", header=header) as modified:

		# 2) Iterate through input SAM/BAM file
		for read in original.fetch(rnam):

			# skip uninteresting alignments
			if (read.is_unmapped) or (read.is_secondary) or (read.is_qcfail): continue

			# declare variables
			#qnam = read.query_name
			qlen = read.query_length
			qseq = read.query_sequence
			lpos = read.reference_start
			rpos = read.reference_end

			# 3) Test if query aligns beyond End of Reference
			if (rpos > rlen) or (lpos < 0): continue #print("{}\t{}\t{}\t{}\t{}".format(rnam,rlen,qnam,lpos,rpos))

			# determine whether read should contain CT or GA
			CT = bool((read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse))

			# get aligned pairs
			alignments = read.get_aligned_pairs()
			mismatches = 0

			# 4) Iterate through each aligned position to identify mismatches
			for bp in alignments:

				# skip indels
				if (bp[0] == None) or (bp[1] == None): continue 

				# count opposite-strand mismatches
				if (qseq[bp[0]] == "A" and rseq[bp[1]] == "G" and CT) or (qseq[bp[0]] == "T" and rseq[bp[1]] == "C" and not CT):
					mismatches += 1
			
			# 5) test if opposite-strand mismatches reach threshold proportion cutoff
			if mismatches > 0:
				
				prop = mismatches/qlen

				if prop > CUTOFF: continue #print("{}\t{}\n{}".format(qnam,mismatches,qseq))
				else: modified.write(read)

			else: modified.write(read)

	# return the filename
	return name



####### Function for merging bam files recursively according to ulimit
def merger(TDIR,batch):

	'''
	TDIR = path to temp dir for temporary merge files eg. "/var/tmp/adfas7d"
	batch = list of files to merge
	'''

	# set the tempfile name
	bam = tempfile.mktemp(dir=TDIR)

	# merge list of files from batch
	arguments = ["-cpf",bam]
	arguments = arguments + batch
	pysam.merge(*arguments)

	# return merged bam
	return bam



####### Function for splitting jobs according to open file limit
def recursive_merger(jobs,rlimit,pool,TDIR):

	# jobs contains a list of files to be merged
	if len(jobs) > 1:
		merges = []

		# get batches of files and submit them to the merger process
		for i in range(0, len(jobs), int(rlimit/2)):
			batch = []
			for job in jobs[i:i+int(rlimit/2)]: batch.append(job.get())
			merge = pool.apply_async(merger, (TDIR, batch))
			merges.append(merge)
		
		# now we have a new list of files so reapply the function
		jobs = recursive_merger(merges,rlimit,pool,TDIR)
	
	# return the final list of jobs when the list contains 1
	return jobs



## END OF FUNCTIONS
###################

#############
## RUN SCRIPT

# define argparse
usage = ''' This program takes an input BAM file, typically from erne-bs5, and filters the
	alignments to exclude beyond-end-of-reference alignments and those that exceed
	a user-defined threshold of opposite-strand mismatches. '''

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('infasta', metavar = 'in.fa', help='[REQUIRED] The path to the reference FASTA file')
parser.add_argument('infile', metavar = 'in.bam', help='[REQUIRED] The path to the original SAM/BAM file')
parser.add_argument('outfile', metavar = 'out.bam', help='[REQUIRED] The path to the output SAM/BAM file')
parser.add_argument('-T','--temp', metavar='<TEMP>', help='[OPTIONAL] Path to temp directory [default: /var/tmp]', default="/var/tmp")
parser.add_argument('-c', '--cutoff', metavar='', help='[OPTIONAL] The threshold proportion of opposite bisulfite mismatches on strand, by which to filter [default: 0.05]', type=float, default=0.05)
parser.add_argument('-t','--threads', metavar='', help='[OPTIONAL] Number of threads for reading [default: 1]', type=int, default=1)

args = parser.parse_args()

# call main()
if __name__ == '__main__':
	main(args.infasta,args.infile,args.outfile,args.temp,args.cutoff,args.threads)

## END OF SCRIPT
################