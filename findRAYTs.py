#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
from micomplete import calcCompleteness
from micomplete import micomplete
from contextlib import contextmanager
from distutils import spawn
from collections import defaultdict
import sys
import re
import os
import argparse
import shutil
import multiprocessing as mp
import subprocess

def worker(fasta, seqType, name, hmm, evalue=1e-20, outfile=sys.stdout):
    if seqType == "faa":
        faa = fasta
    elif seqType == "fna":
        faa = micomplete.create_proteome(fasta, name)
    elif re.match("(gb.?.?)|genbank", seqType):
        faa = micomplete.extract_gbk_trans(fasta, name)
    else:
        raise TypeError('Sequence type needs to be one of faa/fna/gbk')
    baseName = os.path.basename(fasta).split('.')[0]
    if not name:
        name = baseName
    get_RAYTs(faa, name, hmm, evalue)
    return faa

def listener(q):
    """Process to write results in a thread safe manner"""
    handle = open_stdout(outfile)
    while True:
        output = q.get()
        if type(output) == str:
            handle.write()
        else:
            continue

@contextmanager
def open_stdout(outfile):
    """Function dynamically allows printing to file or stdout when used as open 
    with 'with'. No filename or '-' results in stdout, any other name will be
    used as filename to be written"""
    if outfile and outfile != '-':
        handle = open(outfile, 'w')
    else:
        handle = sys.stdout
    # close open file if not stdout
    try:
        yield handle
    finally:
        if handle is not sys.stdout:
            handle.close()

def get_RAYTs(faa, name, hmm, evalue=1e-20):
    """Retrieves markers in hmm from the given proteome"""
    comp = calcCompleteness(faa, name, hmm, evalue)
    foundRAYTs, dupRAYTs, totl = comp.get_completeness()
    # ensure unique, best RAYT match for each gene
    # but allow infinite gene matches for each RAYT
    print(foundRAYTs)
    gene_matches = defaultdict(list)
    for RAYT, match in foundRAYTs.items():
        for gene, evalue in match:
            if gene not in gene_matches:
                gene_matches[gene].append([RAYT, evalue])
            elif float(evalue) < float(gene_matches[gene][0][1]):
                print(evalue)
                gene_matches[gene].pop()
                gene_matches[gene].append([RAYT, evalue])
            print(gene + evalue)
            print(gene_matches[gene])
        print(gene_matches[gene])


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

    with open(args.fastaList) as seq_file:
        inputSeqs = [ seq.strip().split('\t') for seq in seq_file ]
    print(inputSeqs)
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=int(args.threads) + 1)
    # init listener here

    jobs = []
    for i in inputSeqs: 
        print(i)
        if len(i) == 2:
            i.append(None)
        job = pool.apply_async(worker, (i[0], i[1], i[2], args.hmms))
        jobs.append(job)
    # get() all processes to catch errors
    for job in jobs:
        job.get()
    q.put("done")
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
