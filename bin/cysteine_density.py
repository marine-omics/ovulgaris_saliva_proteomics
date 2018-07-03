#!/usr/bin/env python

# Calculating the cysteine density in a moving frame window for FASTA file containing numerous sequences

#to keep the program for general use, argparse and sys used to allow the input pile to be user defined through the command line

import argparse
import sys
import string
#from Bio import SeqIO

#FRAME_ADVANCE = 1;
CYSTEINE = 'C'

parser = argparse.ArgumentParser(description='Calculate cysteine density for protein sequences.')

parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),default=sys.stdin)

#defining options that are able to be set through the command line

parser.add_argument('-m','--min-count', type=int, default=6, help='Minimum number of cysteines within window to call a knot. Default 6')
parser.add_argument('-w','--window', type=int, default=30, help='Size of window over which cysteine counts are calculated')
parser.add_argument('-f','--frame', type=int, default=1, help='Length of moving frame')

args = parser.parse_args()



for record in SeqIO.parse(args.infile, "fasta"):
	seq = str(record.seq)

	knotFound = False;
	max_c_count = 0
#to do: CREATE another loop for more specific cysteine knot as a separate column (add it as e.g + str(superknotFound)) copy below section

	#Looping through sequence, by size FRAME_ADVANCE
	for i in range(0, len(seq), args.frame):
		frame = seq[i:args.window+i]
		cCount = 0;
		#print frame
		cCount = string.count(frame,CYSTEINE)
		#Count how many CYSTEINE are in frame
		# for j in range(len(frame)):
		#     if (CYSTEINE in frame[j]):
		#         cCount = cCount + 1;

		if ( cCount > max_c_count):
			max_c_count = cCount;

		#DOES A KNOT EXIST?
		if (cCount >= args.min_count):
			knotFound = True;

	print record.id + "\t" + str(max_c_count) + "\t" + str(knotFound)

	# if (knotFound):
	#     print "knot found: " + record.id
