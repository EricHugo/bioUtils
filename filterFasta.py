#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
import sys
import re
import argparse

class filterFasta():
    def __init__(self, fasta, headers):
        self.seq = SeqIO.to_dict(SeqIO.parse(fasta, "fasta")) 
        self.seqHeaders = [ header for header, seq in self.seq.items() ]
        self.tarHeaders = [ header.strip() for header in open(headers) ]

    def filter_sequences(self, exclude=False):
        matched = []
        for header in self.tarHeaders:
            r = re.compile(header)
            matched.append(filter(r.match, self.seqHeaders))
        matches = [ ''.join(match) for match in matched ]
        if exclude:
            selected = list(set(self.seqHeaders).difference(set(matches)))
            print("sequences to be excluded: "  + str(len(self.seqHeaders) - 
                len(selected)))
        else:
            selected = matches
            print("sequences to be included: " + str(len(selected)))
        seqObjects = [ self.seq[identifier] for identifier in selected ]
        SeqIO.write(seqObjects, "test", "fasta")
        pass


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="""Filter a fasta sequence based
            on a file of headers, to either only include those sequences or include 
            all sequences except those.""")

    parser.add_argument("fastaFile", help="""Fasta file to be filtered""")
    parser.add_argument("headerFile", help="""List of headers for sequences to 
            either include or exclude""")
    parser.add_argument("-e", "--exclude", default=False, action="store_true", 
            required=False, help="""Flag to exclude given headers rather than
            include""")
    parser.add_argument("-o", "--outfile", required=False, help="""Filename to 
            save resulting sequences to. Otherwise prints to stdout.""")

    args = parser.parse_args()

    f = filterFasta(args.fastaFile, args.headerFile)
    f.filter_sequences(args.exclude)
