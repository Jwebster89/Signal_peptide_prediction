#!/usr/bin/env python3

import os, sys, csv
from Bio import SeqIO

def read_accessions(gff):
	with open(gff) as input_h:
		accessions=[]
		rd = csv.reader(input_h, delimiter="\t", quotechar='"')
		for row in rd:
			accessions.append(row[0])
	return(accessions)

def find_accessions(proteins,acc,outfile):
	with open(outfile, 'w') as out_h:
		for seq_record in SeqIO.parse(proteins, "fasta"):
			if seq_record.id in acc:
				SeqIO.write(seq_record, out_h, 'fasta')

def main():
	acc=read_accessions(file)
	find_accessions(proteins,acc, outfile)

file=sys.argv[1]
proteins=sys.argv[2]
outfile=sys.argv[3]
main()
