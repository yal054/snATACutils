#!/bin/python

import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='filter bedpe based on cluster annotations')
parser.add_argument('--cluster', type=str, dest="cluster", help='input cluster')
parser.add_argument('--indir', type=str, dest="indir", help='directory path of input bedpes')
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
    print("filter out bedpe files")
    generate_bedpes(clusterf, dirPrefix, ncpu, outPrefix)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def generate_bedpes(clusterf, dirPrefix, ncpu, outPrefix):
    clust_dat = pd.read_csv(clusterf, sep="\t")
#    clust_dat.rename(columns={'x.sp@sample':'sample','x.sp@cluster':'cluster',}, inplace=True)
#    clust_dat.rename(columns={'L3cluster':'cluster'}, inplace=True)
    clust_dat["uniq_barcode"] = clust_dat[['sample', 'barcode']].apply(lambda x: '.'.join(x), axis=1)
    clust_dat_dict = pd.Series(clust_dat["cluster"].values,index=clust_dat["uniq_barcode"]).to_dict()
    outf_dict = dict()
#    for cls_id in range(1,np.max(clust_dat["cluster"])+1):
    for cls_id in clust_dat["cluster"]:
        outbedpefname = "".join((outPrefix, ".cluster.", str(cls_id), ".bedpe.gz"))
        outf_dict[cls_id] = outbedpefname
    sample_id = pd.unique(clust_dat["sample"])
    p = Pool(ncpu)
    for sample in sample_id:
        p.apply_async(generate_bedpe_worker, (clust_dat_dict, outf_dict, sample, dirPrefix, outPrefix))
    p.close()
    p.join()

def generate_bedpe_worker(clust_dat_dict, outf_dict, sample, dirPrefix, outPrefix):
    print("extract reads from", sample)
    bedpefname = "".join((dirPrefix, sample, ".bedpe.gz"))
    with gzip.open(bedpefname,'rt') as bedpeF:
        for line in bedpeF:
            l = line.strip().split("\t")
            read_name = line.strip().split("\t")[6]
            old_barcode = read_name.split(":")[0]
            new_barcode = ".".join((sample, old_barcode))
            #wline = line.split("\t")[0] + "\t" + str(line.split("\t")[1]) + "\t" + str(line.split("\t")[2]) + "\t" + new_barcode + "\n"
            if new_barcode in clust_dat_dict:
                #wline = line.split("\t")[0] + "\t" + str(line.split("\t")[1]) + "\t" + str(line.split("\t")[2]) + "\t" + line.split("\t")[3] + "\t" + str(line.split("\t")[4]) + "\t" + str(line.split("\t")[5]) + "\t" + new_barcode + "\t" + line.split("\t")[7] + "\t" + str(line.split("\t")[8]) + "\t" + str(line.split("\t")[9]) + "\n"
                wline = str(l[0]) + "\t" + str(l[1]) + "\t" + str(l[2]) + "\t" + str(l[3]) + "\t" + str(l[4]) + "\t" + str(l[5]) + "\t" + new_barcode + "\t" + str(l[7]) + "\t" + str(l[8]) + "\t" + str(l[9]) + "\n"
                cls = clust_dat_dict[new_barcode]
                outbedpefname = outf_dict[cls]
                obedpe = gzip.open(outbedpefname,'a+')
                obedpe.write(wline.encode())
            else:
                continue
    obedpe.close()
    bedpeF.close()

if __name__ == "__main__":
    """filter bedpe based on cluster annotations"""
    run()


