#!/usr/bin/env python

'''
Title: correct_sam_format.py
Date: 20180730
Author: Adam Nunn
Description:
	Fix the FLAG and CIGAR values in a paired-end SAM file generated by erne-bs5.
	Specifically, correct the custom bisulfite strand information in the FLAG and
	remove any zero-length elements present in the CIGAR string. Input BAM files
    must be sorted and indexed for parallel processing.

    #### Standard ERNE-BS5 FLAG values ####
    #  65 - read paired, first in pair
    # 129 - read paired, second in pair
    # 113 - read paired, read reverse strand, mate reverse strand, first in pair
    # 177 - read paired, read reverse strand, mate reverse strand, second in pair

    #### SAM-spec-compliant FLAG values ####
    #  97 - read paired, mate reverse strand, first in pair
    #  99 - read paired, read mapped in proper pair, mate reverse strand, first in pair
    # 145 - read paired, read reverse strand, second in pair
    # 147 - read paired, read mapped in proper pair, read reverse strand, second in pair
    #  81 - read paired, read reverse strand, first in pair
    #  83 - read paired, read mapped in proper pair, read reverse strand, first in pair
    # 161 - read paired, mate reverse strand, second in pair
    # 163 - read paired, read mapped in proper pair, mate reverse strand, second in pair

List of functions:
	main()
    worker()
    merger()
    recursive_merger()

Procedure:
	1. Declare variables and build threads pool
	2. Open initial BAM file instance to read references
	3. Fire off workers to process reads for each scaffold 'ref' from 'BAM'
	4. Merge results from the workers with a recursive merging function

Usage:
	./correct_sam_format.py [-h, --help] \
        [-i, --insert] int \
        [-t, --threads] int \
        [-T, --temp] <TEMP> \
        [infile] [outfile]
eg. ./correct_sam_format.py -i 500 -t 4 -T /var/tmp in.bam out.bam
'''

###################
## INIT ENVIRONMENT

import multiprocessing as mp
import resource, tempfile, os, argparse
import pysam, re


##################
## DEFINE __MAIN__
def main(BAM,OUT,TEMP,THREADS,MAXINS):

    # 1) Declare variables and build threads pool
    pool = mp.Pool(int(THREADS))
    jobs = []

    # determine open file limit
    rlimit = resource.getrlimit(resource.RLIMIT_NOFILE)
    rlimit = rlimit[0]

    # define temp directory
    TDIR = tempfile.mkdtemp(dir=TEMP)

    # 2) Open initial BAM file instance to read references
    with pysam.AlignmentFile(BAM, "rb") as original:

        header = original.header.to_dict()
        if "CO" in header: del header["CO"]

        for ref in original.get_index_statistics():

            # 3) Fire off workers to process reads for each scaffold 'ref' from 'BAM'
            if ref.mapped > 0:
                job = pool.apply_async(worker, (BAM, TDIR, header, ref.contig, MAXINS))
                jobs.append(job)

    # 4) Merge results from the workers with a recursive merging function
    bam = recursive_merger(jobs,rlimit,pool,TDIR)
    bam = bam[0].get()
    os.replace(bam, OUT)

    # close the pool
    pool.close()
    pool.join()


###################
## DEFINE FUNCTIONS

####### Function for READING the input reads, modifying, and sending them to 'q'
def worker(BAM,TDIR,header,rnam,MAXINS):

    '''
    BAM = path to input bam file eg. "/path/to/input.bam"
    TDIR = path to temp dir for processed bam files eg. "/var/tmp/adfas7d"
    header = dictionary object containing modified header for writing
    ref = current scaffold or chromosome reference eg. "Chr1"
    MAXINS = maximum insert size integer eg. 500
    '''

    # generate a name for the modified bam file
    name = TDIR + "/" + rnam + ".bam"

    # Open pysam.AlignmentFile objects for reading and writing
    with pysam.AlignmentFile(BAM, "rb") as original, pysam.AlignmentFile(name, "wb", header=header) as modified:

        # declare initial vars
        mates = dict()

        # Iterate through input SAM/BAM file
        for read in original.fetch(rnam):

            # skip uninteresting alignments
            if (read.is_unmapped == True) or (read.is_secondary == True) or (read.is_qcfail == True): continue

            # replace zero-length elements
            if read.cigarstring: read.cigarstring = re.sub("\D0","",read.cigarstring)

            # add similar tags to segemehl
            HI = read.get_tag('HI') - 1
            read.set_tag('HI',HI,value_type='i')
            read.set_tag('XB','XX/XX',value_type='Z')

            # the 'read' and 'mate' alignments should not overwrite each other here
            # both pair1pos and pair2pos will be set accordingly when pair is written to file
            if read.is_read1: pair1pos = read.reference_end
            elif read.is_read2: pair2pos = read.reference_end
            
            # set qnam for lookup
            qnam = read.query_name

            ### current read is First of a new pair, so save it as 'mate' and move on to the next
            if not (qnam in mates):
                fragment = MAXINS + read.query_length
                if read.is_read1: mates[qnam] = (read, pair1pos, fragment)
                else: mates[qnam] = (read, pair2pos, fragment)
                continue
			
            ### current read has a mate!!
            else:

                ### pop 'mate' from the dictionary
                if read.is_read1: mate, pair2pos, fragment = mates.pop(qnam)
                else: mate, pair1pos, fragment = mates.pop(qnam)

                ### set fragment size
                fragment = fragment + read.query_length

                ### set flags for 'read' + 'mate'
                for pair in [mate,read]:

                    # read paired, FIRST in pair ->
                    # read paired, FIRST in pair, mate REVERSE strand
                    if pair.flag == 65:
                        if (pair.reference_start < pair.next_reference_start) and \
                            ((pair2pos - pair.reference_start) < fragment): pair.flag = 99
                        else: pair.flag = 97

                    # read paired, SECOND in pair ->
                    # read paired, SECOND in pair, read REVERSE strand
                    elif pair.flag == 129:
                        if (pair.next_reference_start < pair.reference_start) and \
                            ((pair.reference_end - pair.next_reference_start) < fragment): pair.flag = 147
                        else: pair.flag = 145

                    # read paired, read REVERSE strand, mate REVERSE strand, FIRST in pair ->
                    # read paired, read REVERSE strand, FIRST in pair
                    elif pair.flag == 113:
                        if (pair.reference_end > pair2pos) and \
                            ((pair.reference_end - pair.next_reference_start) < fragment): pair.flag = 83
                        else: pair.flag = 81

                    # read paired, read REVERSE strand, mate REVERSE strand, SECOND in pair ->
                    # read paired, mate REVERSE strand, SECOND in pair
                    elif pair.flag == 177:
                        if (pair1pos > pair.reference_end) and \
                            ((pair1pos - pair.reference_start) < fragment): pair.flag = 163
                        else: pair.flag = 161

                    # we're not interested in miscellaneous flag settings for PE reads
                    else: continue

                    # write current 'pair' to output file
                    modified.write(pair)

        # now that file has been iterated, remaining reads in 'mates' are singletons
        for mate in mates:

            # get the alignment and write unchanged to file
            read = mates[mate][0]
            modified.write(read)

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

	# merge
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
usage = ''' Fix the FLAG and CIGAR values in a paired-end SAM file generated by erne-bs5.
	Specifically, correct the custom bisulfite strand information in the FLAG and 
	remove any zero-length elements present in the CIGAR string. '''

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('infile', metavar = 'in.bam', help='[REQUIRED] The path to the old SAM/BAM file')
parser.add_argument('outfile', metavar = 'out.bam', help='[REQUIRED] The path to the new SAM/BAM file')
parser.add_argument('-T','--temp', metavar='<TEMP>', help='[OPTIONAL] Path to temp directory [default: /var/tmp]', default="/var/tmp")
parser.add_argument('-t','--threads', metavar='', help='[OPTIONAL] Number of threads for reading [default: 1]', type=int, default=1)
parser.add_argument('-i','--insert', metavar='', help='[OPTIONAL] The maximum insert size for PE reads [default: 500]', type=int, default=500)

args = parser.parse_args()

# call main()
if __name__ == '__main__':
	main(args.infile,args.outfile,args.temp,args.threads,args.insert)

## END OF SCRIPT
################