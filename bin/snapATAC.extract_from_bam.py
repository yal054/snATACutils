#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='filter bam based on QNAMES')
parser.add_argument('--cluster', type=str, dest="cluster", help='input cluster')
parser.add_argument('--indir', type=str, dest="indir", help='directory path of input bam')
parser.add_argument('--outprfx', type=str, dest="outprfx", help='output directory')

args = parser.parse_args()

import numpy as np
import pandas as pd
import pysam
from time import perf_counter as pc

def run():
    """ Run standard NMF on rank """
    start_time = pc()
    """ init input files """
    clusterf = args.cluster
    dirPrefix = str(args.indir)
    outPrefix = args.outprfx
    print("filter out bam files")
    generate_bams(clusterf, dirPrefix, outPrefix)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def generate_bams(clusterf, dirPrefix, outPrefix):
    clust_dat = pd.read_csv(clusterf, sep="\t")
    clust_dat["uniq_barcode"] = clust_dat[['sample', 'barcode']].apply(lambda x: '.'.join(x), axis=1)
    clust_dat_dict = pd.Series(clust_dat["cluster"].values,index=clust_dat["uniq_barcode"]).to_dict()
    outf_dict = dict()
    sample_id = pd.unique(clust_dat["sample"])
    templateBamF = pysam.AlignmentFile("".join((dirPrefix, sample_id[0], ".bam")))
    for cls_id in pd.unique(clust_dat["cluster"]):
        outsamfname = "".join((outPrefix, ".cluster.", str(cls_id), ".sam"))
        opysam = pysam.AlignmentFile(outsamfname, "wh", template=templateBamF)
        outf_dict[cls_id] = opysam
    if len(sample_id) == 1:
        sample = sample_id[0]
        generate_bam_worker(clust_dat_dict, outf_dict, sample, dirPrefix, outPrefix)
    else:
        for sample in sample_id:
             generate_bam_worker(clust_dat_dict, outf_dict, sample, dirPrefix, outPrefix)

def generate_bam_worker(clust_dat_dict, outf_dict, sample, dirPrefix, outPrefix):
    print("extract reads from", sample)
    bamfname = "".join((dirPrefix, sample, ".bam"))
    bamF = pysam.AlignmentFile(bamfname)
    for b in bamF.fetch(until_eof=True):
        old_barcode = b.query_name.split(':')[0]
        new_barcode = ".".join((sample, old_barcode))
        if new_barcode in clust_dat_dict:
            rest_qname = ":".join(b.query_name.split(':')[1:])
            new_qname = ":".join((str(new_barcode),rest_qname))
            b.query_name = str(new_qname)
            cls = clust_dat_dict[new_barcode]
            osam = outf_dict[cls]
            osam.write(b)
        else:
            continue
    osam.close()
    bamF.close()

if __name__ == "__main__":
    """filter bam based on QNAMES"""
    run()


