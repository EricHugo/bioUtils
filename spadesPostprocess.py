#!/usr/bin/env python

from Bio import SeqIO
import sys

if __name__ == "__main__":
    assembly = sys.argv[1]
    try:
        cutoff = float(sys.argv[2])
    except IndexError:
        cutoff = 500
    try:
        coverage = float(sys.argv[3])
    except IndexError:
        coverage = 10
    p_assembly = SeqIO.parse(assembly, 'fasta')
    processed = []
    for each in p_assembly:
        cov = float(each.name.split('_')[-1])
        if cov < coverage:
            continue
        length = float(each.name.split('_')[-3])
        if length < cutoff:
            continue
        processed.append(each)
        SeqIO.write(each, sys.stdout, 'fasta')
