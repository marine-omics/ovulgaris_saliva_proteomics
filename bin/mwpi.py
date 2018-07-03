#!/usr/bin/env python

# Calculating the molecular weight and isoelectric point for a FASTA file containing numerous sequences

#to keep the program for general use, argparse and sys used to allow the input pile to be user defined through the command line

import argparse
import sys

parser = argparse.ArgumentParser(description='Calculate mw and pi for protein sequences.')

parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),default=sys.stdin)
args = parser.parse_args()

#to read from a FASTA file with a loop over entries using SeqIO define the FASTA sequences and analyse them by ProteinAnalysis
#display the sequence names, molecular weight and isoelectric point

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
for record in SeqIO.parse(args.infile, "fasta"):
	seq = str(record.seq)
	try:
		my_prot = ProteinAnalysis(seq)
		print '{}\t {}\t {}'.format(record.id, my_prot.molecular_weight(), my_prot.isoelectric_point())
	except:
		print '{}\t NaN\t NaN'.format(record.id)
