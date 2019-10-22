#!/usr/bin/env python

'''
Title: filter_sam_uniqs.py
Date: 20180730
Author: Adam Nunn
Description:
  Read a SAM/BAM file and filter multi-mappings based on NH flag, directing uniquely mapping reads to
  one output file and multi-mapping reads to another output file. 

List of functions:
  main()

Procedure:
  1. Open pysam.AlignmentFile objects for reading and writing
  2. Iterate through input SAM/BAM file
  3. Test alignment NH tag for values of 1

Usage:
    ./filter_sam_uniqs.py [infile] [uniqfile] [multfile]
eg. ./filter_sam_uniqs.py in.bam uniq.bam mult.bam
'''

###################
## INIT ENVIRONMENT

import argparse
import pysam

##################
## DEFINE __MAIN__
def main(BAM,UNIQ,MULT):

	# 1) Open pysam.AlignmentFile objects for reading and writing
	with pysam.AlignmentFile(BAM, "rb") as original, pysam.AlignmentFile(UNIQ, "wb", header=original.header) as uniq, \
		pysam.AlignmentFile(MULT, "wb", header=original.header) as mult:

		# 2) Iterate through input SAM/BAM file
		for read in original:

			# 3) Test alignment NH tag for values of 1
			NH = read.get_tag('NH')
			if NH == 1: uniq.write(read)
			else: mult.write(read)


## END OF __MAIN__
##################

#############
## RUN SCRIPT

# define argparse
usage = 'Read a SAM/BAM file and filter multi-mappings based on NH flag, directing uniquely mapping reads to one output file and multi-mapping reads to another output file.'

parser = argparse.ArgumentParser(description=usage)

parser.add_argument('infile', metavar = 'in.bam', help='[REQUIRED] The path to the original SAM/BAM file')
parser.add_argument('ufile', metavar = 'uniq.bam', help='[REQUIRED] The path to the output SAM/BAM file containing unique mapping alignments')
parser.add_argument('mfile', metavar = 'mult.bam', help='[REQUIRED] The path to the output SAM/BAM file containing multi-mapping alignments')

args = parser.parse_args()

# call main()
if __name__ == '__main__':
	main(args.infile,args.ufile,args.mfile)

## END OF SCRIPT
################