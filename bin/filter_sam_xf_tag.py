#!/usr/bin/env python

'''
Title: filter_sam_xf_tag.py
Date: 20180918
Author: Adam Nunn
Description:
	Take an original SAM/BAM file and filter based on a user-defined proportion
	of opposite bisulfite mismatches on each strand, indicating possible genomic
	inversions in the query sample. Reads containing XF tags with a higher
	proportion than the threshold will be excluded from the output SAM/BAM file
	[default: 0.05].

List of functions:
  main()

Procedure:
  1. Open pysam.AlignmentFile objects for reading and writing
  2. Iterate through input SAM/BAM file
  3. Test for presence of XF tag
  4. Skip alignment if XF proportion is greater than user-defined threshold

Usage:
    ./filter_sam_xf_tag.py [-c, --cutoff] <float> [infile] [outfile]
eg. ./filter_sam_xf_tag.py -c 0.05 in.bam out.bam
'''

###################
## INIT ENVIRONMENT

import argparse
import pysam

##################
## DEFINE __MAIN__
def main(BAM,OUT,CUTOFF):

	# 1) Open pysam.AlignmentFile objects for reading and writing
	with pysam.AlignmentFile(BAM, "rb") as original, pysam.AlignmentFile(OUT, "wb", header=original.header) as modified:

		# 2) Iterate through input SAM/BAM file
		for read in original:

			# skip uninteresting reads
			if (read.is_unmapped) or (read.is_secondary) or (read.is_qcfail): continue
	
			# 3) Test for presence of XF tag
			if read.has_tag('XF'):

				XF = read.get_tag('XF')
				qlen = read.query_length
			
				# Test for non-zero 
				if XF > 0:

					# 4) Skip alignment if XF proportion is greater than user-defined threshold
					prop = XF/qlen
					if prop > CUTOFF: continue #print("{}\t{}\t{}\t{}".format(read.query_name,strand,XF,seq))

				# write read to file
				modified.write(read)

			# no XF tags present in read, so quit
			else:
				print("ERROR: XF tags are missing from the input SAM/BAM file")
				raise SystemExit(1)


## END OF __MAIN__
##################

#############
## RUN SCRIPT

# define argparse
usage = 'Filter SAM/BAM file based on user-defined proportion of opposite-strand mismatches.'

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('infile', metavar = 'in.bam', help='[REQUIRED] The path to the original SAM/BAM file')
parser.add_argument('outfile', metavar = 'out.bam', help='[REQUIRED] The path to the modified SAM/BAM file')
parser.add_argument('-c', '--cutoff', metavar='', help='[OPTIONAL] The threshold proportion of opposite bisulfite mismatches on strand, by which to filter [default: 0.05]', type=float, default=0.05)

args = parser.parse_args()

# call main()
if __name__ == '__main__':
	main(args.infile,args.outfile,args.cutoff)

## END OF SCRIPT
################