#!/usr/bin/env python

# Countthe amino acids for a FASTA file containing numerous sequences

#to keep the program for general use, argparse and sys used to allow the input pile to be user defined through the command line

import argparse
import sys

parser = argparse.ArgumentParser(description='Calculate mw and pi for protein sequences.')

parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),default=sys.stdin)
args = parser.parse_args()

#to read from a FASTA file with a loop over entries using SeqIO define the FASTA sequences and analyse them by ProteinAnalysis
#display the sequence names, molecular weight and isoelectric point


from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
for record in SeqIO.parse(args.infile, "fasta"):
    seq = str(record.seq)
#    my_c = Seq(seq)
    my_prot = ProteinAnalysis(seq)
    aa_counts = my_prot.count_amino_acids()
    c_counts = aa_counts['C']
    print '{}\t {}'.format(record.id,c_counts)

    #print '{}\t {}'.format(record.id, my_c.count("C"))
