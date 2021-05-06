#!/usr/bin/env

import sys
from Bio import SeqIO


def find_dupes(seq, seq_list):
    seq_headers = [seqi.id for seqi in SeqIO.parse(seq, "genbank")]
    for each in seq_list:
        headers = [seqi.id for seqi in SeqIO.parse(each, "genbank")]
        for header in seq_headers:
            if header in headers:
                print("{} match in {}".format(seq, each))


if __name__ == "__main__":
    seq_f = sys.argv[1]
    seq_list = []
    with open(seq_f) as f:
        for seq in f:
            seq_list.append(seq.strip())
    while seq_list:
        seq = seq_list.pop()
        find_dupes(seq, seq_list)

