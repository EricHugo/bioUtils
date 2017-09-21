#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
from micomplete import calcCompleteness
from micomplete import micomplete
import sys
import re
import argparse
import shutil
import spawn
import multiprocessing as mp
import subprocess

def worker(fasta,  name, hmm, evalue=1e-20):
    pass

def get_RAYTs(fasta, name, hmm, evalue=1e-20):
    if not fasta == "faa":
        faa = micomplete.extract_gbk_trans(fasta, name)
        faa = micomplete.create_proteome(fasta, name)
    else:
        faa = fasta
    comp = calcCompleteness(fasta, name, hmm, evalue)
    
    pass


def main():
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
    parser.add_argument("--threads", required=False, default=1, 
            help="""Number of threads to be used concurrently""")

    args = parser.parse_args()

    try:
        assert shutil.which('hmmsearch')
    except AssertionError:
        raise RuntimeError('Unable to find hmmsearch in path')
    except AttributeError:
        try:
            assert spawn.find_executable('hmmsearch')
        except AssertionError:
            raise RuntimeError('Unable to find hmmsearch in path')

    with open(args.sequence) as seq_file:
        inputSeqs = [ seq.strip().split('\t') for seq in seq_file ]

    manager = mp.Manger()
    q = manager.Queue()
    pool = mp.Pool(process=args.threads + 1)
    # init listener here

    jobs = []
    for i in inputSeqs: 
        if len(i) == 2:
            i.append(None)
        job = pool.apply_async(worker(i[0], i[1], i[2]))
        jobs.append(job)
    # get() all processes to catch errors
    for job in jobs:
        job.get()
    q.put("done")
    weights_file = writer.get()
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
