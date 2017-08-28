#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
import Bio
import sys
import argparse

class translateFasta():
    def __init__(self, fasta): 
        self.seq = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    def translate_sequences(self, outfile=None, ttable=11):
        if outfile:
            output_handle = open(outfile, mode='w+')
        else:
            output_handle = sys.stdout
        for identifier, seqObject in self.seq.items():
            output_handle.write(">" + identifier + '\n')
            try:
                output_handle.write(str(seqObject.seq.translate(to_stop=True,
                    table=ttable, cds=True)) + '\n')
            except Bio.Data.CodonTable.TranslationError:
                print("Warning: sequence %s not considered CDS. Writing \
                        translation regardless." % identifier, file=sys.stderr)
                output_handle.write(str(seqObject.seq.translate(to_stop=True,
                    table=ttable)) + '\n')


if __name__ == "__main__":
        
    parser = argparse.ArgumentParser(description="""Simple translation of 
            given fasta file. Warns if not valid CDS""")

    parser.add_argument("fastaFile", help="""Fasta file to be filtered""")
    parser.add_argument("-o", "--outfile", required=False, help="""Filename to 
                save resulting translations to. Otherwise prints to stdout.""")
    parser.add_argument("-t", "--ttable", required=False, default=11,
                help="""Provide translation table number""")

    args = parser.parse_args()

    translate = translateFasta(args.fastaFile)
    translate.translate_sequences(outfile=args.outfile, ttable=args.ttable)
