#!/bin/python

import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='filter bed based on cluster annotations')
parser.add_argument('--cluster', type=str, dest="cluster", help='input cluster')
parser.add_argument('--indir', type=str, dest="indir", help='directory path of input beds')
parser.add_argument('--cpu', type=int, dest="cpu", default=1, help='# of cpu')
parser.add_argument('--outprefix', type=str, dest="outprefix", help='output prefix')

args = parser.parse_args()

import numpy as np
import pandas as pd
import pysam
from time import perf_counter as pc
import gzip

def run():
    """ Run """
    start_time = pc()
    """ init input files """
    clusterf = args.cluster
    dirPrefix = str(args.indir)
    ncpu = args.cpu
    outPrefix = args.outprefix
    print("filter out bed files")
    generate_beds(clusterf, dirPrefix, ncpu, outPrefix)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def generate_beds(clusterf, dirPrefix, ncpu, outPrefix):
    clust_dat = pd.read_csv(clusterf, sep="\t")
    clust_dat.rename(columns={'x.sp@sample':'sample','x.sp@cluster':'cluster',}, inplace=True)
    clust_dat["uniq_barcode"] = clust_dat[['sample', 'barcode']].apply(lambda x: '.'.join(x), axis=1)
    clust_dat_dict = pd.Series(clust_dat["cluster"].values,index=clust_dat["uniq_barcode"]).to_dict()
    outf_dict = dict()
#    for cls_id in range(1,np.max(clust_dat["cluster"])+1):
    for cls_id in clust_dat["cluster"]:
        outbedfname = "".join((outPrefix, ".cluster.", str(cls_id), ".bed.gz"))
        outf_dict[cls_id] = outbedfname
    sample_id = pd.unique(clust_dat["sample"])
    p = Pool(ncpu)
    for sample in sample_id:
        p.apply_async(generate_bed_worker, (clust_dat_dict, outf_dict, sample, dirPrefix, outPrefix))
    p.close()
    p.join()

def generate_bed_worker(clust_dat_dict, outf_dict, sample, dirPrefix, outPrefix):
    print("extract reads from", sample)
    bedfname = "".join((dirPrefix, sample, ".snap.bed.gz"))
    with gzip.open(bedfname,'rt') as bedF:
        for line in bedF:
            old_barcode = line.strip().split("\t")[3]
            new_barcode = ".".join((sample, old_barcode))
            wline = line.split("\t")[0] + "\t" + str(line.split("\t")[1]) + "\t" + str(line.split("\t")[2]) + "\t" + new_barcode + "\n"
            if new_barcode in clust_dat_dict:
                cls = clust_dat_dict[new_barcode]
                outbedfname = outf_dict[cls]
                obed = gzip.open(outbedfname,'a+')
                obed.write(wline.encode())
            else:
                continue
    obed.close()
    bedF.close()

if __name__ == "__main__":
    """filter bed based on cluster annotations"""
    run()


