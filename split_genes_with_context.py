#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import subprocess
import re
import multiprocessing as mp
from textwrap import wrap
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from termcolor import cprint

def find_gene_seq(genbank_file, locid, upstream, downstream, feature="locus_tag"):
    #hack to fix misc_RNA features having no "old_locus_tag"
    print(locid)
    found = False
    feature_type = "CDS"
    if re.search("misc", locid):
        feature = "locus_tag"
        feature_type = "misc_RNA"
    for seq_record in SeqIO.parse(args.genbank_file, 'genbank'):
        for seq_feature in seq_record.features:
            if seq_feature.type == feature_type and seq_feature.qualifiers[feature][0] == locid:
                found = True
                try:
                    assert seq_feature.qualifiers['gene']                                                             
                except KeyError:                    
                    seq_feature.qualifiers['gene'] = ["no_gene_name"]                                                 
                    #continue                       
                #print(seq_record.seq[seq_feature.location.nofuzzy_start:seq_feature.location.nofuzzy_end])
                fwd = seq_feature.location.nofuzzy_start - upstream
                if fwd < 1:
                    #TODO: allow wrapping back to end of seq
                    print("fwd: " + locid)
                    fwd = 1
                rev = seq_feature.location.nofuzzy_end + downstream
                if rev > len(seq_record.seq):
                    print("rev: " + locid)
                    #TODO: allow wrapping to start of seq
                    rev = len(seq_record.seq)
                # here loop starting positions find any < than bp to preseve upstream
                # gather features from there
                # same for downstream
                # extract seq from upstream-start:downstream-start
                # loop to edit their positions <pos> - (<upstream-start> - 1) (?) 
                features, up, down = find_gene_by_pos(genbank_file, locid, fwd, rev,
                                            seq_feature.location.nofuzzy_start,
                                            seq_feature.location.nofuzzy_end)
                if rev > down:
                    down = rev
                if fwd < up:
                    up = fwd
                gene = seq_feature.qualifiers['gene'][0]
                name = seq_record.name
                seq = seq_record.seq[up:down]
                q.put((features, locid, seq))
    if not found:
        cprint("Could not find: " + locid, "red")
    return 


def find_gene_by_pos(seq, gene, up_lim, down_lim, up_loc, down_loc):
    features = []
    up_locs = []
    down_locs = []
    for seq_record in SeqIO.parse(args.genbank_file, 'genbank'):
        for seq_feature in seq_record.features:
            if seq_feature.location.nofuzzy_end < down_lim and seq_feature.location.nofuzzy_start > up_lim:
                up_locs.append(seq_feature.location.nofuzzy_start)
                down_locs.append(seq_feature.location.nofuzzy_end)
                features.append(seq_feature)
                continue
            if seq_feature.location.nofuzzy_start < down_lim and seq_feature.location.nofuzzy_start > up_lim:
                up_locs.append(seq_feature.location.nofuzzy_start)
                down_locs.append(seq_feature.location.nofuzzy_end)
                features.append(seq_feature)
                continue
            if seq_feature.location.nofuzzy_end > up_lim and seq_feature.location.nofuzzy_end < down_lim:
                up_locs.append(seq_feature.location.nofuzzy_start)
                down_locs.append(seq_feature.location.nofuzzy_end)
                features.append(seq_feature)
                continue
    # incase gene was found to start before up_lim, set that to be the shift_factor
    if up_lim < min(up_locs):
        shift_factor = up_lim
    else:
        shift_factor = min(up_locs)
    shifted_features = []
    # shift to match new sequence location
    for feature in features:
        feature = feature._shift(-shift_factor)
        shifted_features.append(feature)
    return shifted_features, min(up_locs), max(down_locs)


def _listener(q, sid, name, description, dxrefs):
    all_fna = []
    while True:
        features, locid, seq = q.get()
        if locid == "_":
            break
        fname = locid + '.gbk'
        print(fname)
        f = SeqRecord(seq, sid, name, description, dxrefs, features)
        SeqIO.write(f, fname, "genbank")
        continue


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Seperates Genbank file genes\
            from csv list with desired number of bps up- and downstream")

    parser.add_argument("genbank_file", help="Genbank file to be seperated")
    parser.add_argument("gene_csv", help="csv containing list of genes to seperate")
    parser.add_argument("-u", "--upstream", type=int, default=0, help="Number of "\
                        "bp to preserve upstream of gene")
    parser.add_argument("-d", "--downstream", type=int, default=0, help="Number of "\
                        "bp to preserve downstream of gene")
    parser.add_argument("--old_locus", default=False, action='store_true', 
                        help="Flag: use old locus tag for gene matching instead. "\
                        "Default=False")
    parser.add_argument("--threads", required=False, default=1, type=int,
                        help="""Specify number of threads to be used in
                        parallel. Default = 1""")

    args = parser.parse_args()
    
    locids = []
    with open(args.gene_csv) as f:
        for locid in f:
            locids.append(locid.split(',')[0].strip())
    #df = pd.read_csv(args.gene_csv)
    #sys.exit()
    all_fna = []
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads + 1)
    for seq_record in SeqIO.parse(args.genbank_file, 'genbank'):
        writer = pool.apply_async(_listener, (q, seq_record.id, seq_record.name, seq_record.description, seq_record.dbxrefs))
        break
    jobs = []
    for locid in locids:
        job = pool.apply_async(find_gene_seq, (args.genbank_file, locid, args.upstream,
                                               args.downstream))
        jobs.append(job)
    for job in jobs:
        job.get()
    q.put(('_','_','_'))
    pool.close()
    pool.join()
