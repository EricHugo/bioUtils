#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
from micomplete import calcCompleteness
from micomplete import micomplete
import sys
import re
import argparse

def get_RAYTs(fasta, name, hmm, evalue=1e-20):
    comp = calcCompleteness(fasta, name, hmm, evalue)

    pass


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="""For a given set of genomes
            investigate presence of RAYT via HMMs. Attempt to categorize mathces
            and investigate flanking genes and 16mers.""")

    parser.add_argument("fastaList", help="""Sequence(s) along with type (fna, 
                faa, gbk) provided in a tabular format""")
    parser.add_argument("-o", "--outfile", required=False, help="""Filename to 
            save results. Otherwise prints to stdout.""")
    parser.add_argument("--hmms", required=False, default=False,
                help="""Specifies a set of HMMs to be used for completeness check 
                        or linkage analysis""")

    args = parser.parse_args()

