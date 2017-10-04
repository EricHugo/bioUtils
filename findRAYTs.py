#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
from micomplete import calcCompleteness
from micomplete import micomplete
from contextlib import contextmanager
from distutils import spawn
from collections import defaultdict
from ftplib import FTP
from datetime import datetime
import sys
import mmap
import re
import os
import argparse
import shutil
import multiprocessing as mp
import subprocess
import tarfile
import hashlib

class parse_taxonomy():
    def __init__(self, tax_path='taxonomy'):
        self.tax_path = tax_path
        self.download_taxonomy()
        self.tax_names = tax_path + '/names.dmp'
        self.tax_nodes = tax_path + '/nodes.dmp'

    def download_taxonomy(self):
        """Module downloads, unpacks, and verifies ncbi taxdump"""
        try:
            modified_time = datetime.fromtimestamp(os.path.getmtime(self.tax_path +
                '/nodes.dmp'))
            age = datetime.today() - modified_time
            if not age.days >= 2:
                return
        except FileNotFoundError:
            pass
        with FTP('ftp.ncbi.nih.gov') as ncbi:
            ncbi.login()
            ncbi.cwd('pub/taxonomy/')
            print("Downloading taxdump from NCBI...", file=sys.stderr, end='',
                    flush=True)
            ncbi.retrbinary('RETR taxdump.tar.gz', open('taxdump.tar.gz',
                'wb').write)
            ncbi.retrbinary('RETR taxdump.tar.gz.md5', open('taxdump.tar.gz.md5',
                'wb').write)
            print(" Done", file=sys.stderr)
            ncbi.quit()
        try:
            os.mkdir(self.tax_path)
        except FileExistsError:
            pass
        local_md5 = hashlib.md5(open('taxdump.tar.gz', 'rb').read()).hexdigest()
        with open('taxdump.tar.gz.md5') as md5:
            ncbi_md5 = md5.readline().split()[0]
            if not ncbi_md5 == local_md5:
                raise ValueError('Local md5 checksum does not equal value from ncbi')
        tar = tarfile.open('taxdump.tar.gz')
        tar.extractall(path=self.tax_path)
        os.remove('taxdump.tar.gz')
        os.remove('taxdump.tar.gz.md5')

    def find_taxid(self, name):
        """Given full scientific name: return taxid"""
        with open(self.tax_names, 'r') as names_file:
            for tax in names_file:
                if re.fullmatch(name, tax.split('|')[1].strip()):
                    taxid = tax.split('|')[0].strip()
        return taxid

    def find_scientific_name(self, taxid):
        """Given taxid, return full scientific name"""
        with open(self.tax_names, 'r') as names_file:
            for tax in names_file:
                if re.fullmatch(taxid, tax.split('|')[0].strip()):
                    name = tax.split('|')[1]
        return name

    def parse_taxa(self, taxid):
        """Will locate parent of given taxid, loop to find full taxonomic list,
        returns None when super kingdom is reached"""
        print(taxid)
        taxids = taxid.split()
        ranks = []
        with open(self.tax_nodes, 'r+b') as node_file:
            m = mmap.mmap(node_file.fileno(), 0, prot=mmap.PROT_READ)
            while True:
                # dict comp to quickly find match first field and retrieve second
                parent = [ (taxids.append(str(node, "utf-8").split('|')[1].strip()),
                        ranks.append(str(node, "utf-8").split('|')[2].strip()))
                        for node in iter(m.readline, bytes()) 
                        if re.fullmatch(str(taxids[-1]), str(node, "utf-8").split('|')[0].strip()) ]
                print(taxids)
                print(ranks)
                m.seek(0)
                if taxids[-1] == '1':
                    taxids.pop()
                    ranks[0] = "Name"
                    break
        print(taxid)
        return parent
                

def _worker(fasta, seqType, name, hmm, evalue=1e-20, outfile=sys.stdout):
    tax = parse_taxonomy()
    if seqType == "faa":
        faa = fasta
    elif seqType == "fna":
        faa = micomplete.create_proteome(fasta)
    elif re.match("(gb.?.?)|genbank", seqType):
        name = get_gbk_feature(fasta, 'organism')
        faa = micomplete.extract_gbk_trans(fasta, name + '.faa')
        taxid = tax.find_taxid(name)
        tax.parse_taxa(taxid)
    else:
        raise TypeError('Sequence type needs to be one of faa/fna/gbk')
    baseName = os.path.basename(fasta).split('.')[0]
    if not name:
        name = baseName
    gene_matches = get_RAYTs(faa, name, hmm, evalue)
    rayt_gene_list = extract_hits(faa, seqType, gene_matches)
    return faa

def listener(q):
    """Process to write results in a thread safe manner"""
    handle = open_stdout(outfile)
    while True:
        output = q.get()
        if type(output) == str:
            handle.write(output)
        elif output == "done":
            break
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
    # close and flush open file if not stdout
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
    print(gene_matches)
    return gene_matches

def extract_hits(faa, seqType, gene_RAYT):
    """Extract the protein sequences from a given proteome given a dict of 
    gene names"""
    rayt_gene_list = set( gene for gene, rayt in gene_RAYT.items() )
    rayt_proteins = [ seq for seq in SeqIO.parse(faa, "fasta") if seq.id in
            rayt_gene_list ]
    print(rayt_proteins)
    return rayt_proteins

def get_gbk_feature(handle, feature_type):
    """Get specified organism feature from gbk file"""
    input_handle = open(handle, mode='r')
    for record in SeqIO.parse(input_handle, "genbank"):
        for feature in record.features:
            if feature.type == "source":
                value = ''.join(feature.qualifiers[feature_type])
    return value

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

    parse_taxonomy()

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
        job = pool.apply_async(_worker, (i[0], i[1], i[2], args.hmms))
        jobs.append(job)
    # get() all processes to catch errors
    for job in jobs:
        job.get()
    q.put("done")
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
