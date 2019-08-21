#!/bin/python

import argparse
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='filter bed based on cluster annotations')
parser.add_argument('--cluster', type=str, dest="cluster", help='input cluster')
parser.add_argument('--inbed', type=str, dest="inbed", help='input bed')
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
	inbed = str(args.inbed)
	outPrefix = args.outprefix
	print("filter out bed files")

	clust_dat = pd.read_csv(clusterf, sep="\t")
	clust_dat.rename(columns={'x.sp@sample':'sample','x.sp@cluster':'cluster',}, inplace=True)
	clust_dat["uniq_barcode"] = clust_dat[['sample', 'barcode']].apply(lambda x: '.'.join(x), axis=1)
	clust_dat_dict = pd.Series(clust_dat["cluster"].values,index=clust_dat["uniq_barcode"]).to_dict()

	outf_dict = dict()
	for cls_id in range(1,np.max(clust_dat["cluster"])+1):
		outbedfname = "".join((outPrefix, ".cluster.", str(cls_id), ".bed.gz"))
		outf_dict[cls_id] = outbedfname

	with gzip.open(inbed,'rt') as bedF:
		for line in bedF:
			barcode = line.strip().split("\t")[3]
			if barcode in clust_dat_dict:
				cls = clust_dat_dict[barcode]
				outbedfname = outf_dict[cls]
				obed = gzip.open(outbedfname,'a+')
				obed.write(line.encode())
			else:
				continue
	obed.close()
	bedF.close()

	end_time = pc()
	print('Used (secs): ', end_time - start_time)


if __name__ == "__main__":
	"""filter bed based on cluster annotations"""
	run()


