#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='convert 3-cols to sparse matrix')
parser.add_argument('-i', '--in', type=str, dest="input", help='input 3-col file: peak, cell, ct')
parser.add_argument('-o', '--out', type=str, dest="output", help='prefix of output file')

args = parser.parse_args()

import numpy as np
from scipy import sparse
from scipy.sparse import save_npz
from time import perf_counter as pc

def run():
	""" Run and count to peak """
	start_time = pc()
	""" init input files """
	inf = args.input
	outf = args.output
	cells = {}
	peaks = {}
	with open(inf) as in_col3dat:
		xgi_list = []
		ygi_list = []
		ct_list = []
		ncell = 0
		npeak = 0
		for idx,line in enumerate(in_col3dat):
			peak = line.split()[0]
			cell = line.split()[1]
			ct = line.split()[2]
			if peak not in peaks:
				peaks[peak] = npeak
				npeak += 1
			if cell not in cells:
				cells[cell] = ncell
				ncell += 1
			xgi_list.append(cells[cell])
			ygi_list.append(peaks[peak])
			ct_list.append(int(ct))
	in_col3dat.close()
	cooMx = sparse.coo_matrix((ct_list, (xgi_list,ygi_list)))
	outxgi = ".".join([outf,"xgi"])
	outygi = ".".join([outf,"ygi"])
	outnpz = ".".join([outf,"npz"])
	save_npz(outnpz, cooMx)
	with open(outxgi, "w") as xgif:
		for key in cells:
			print(key, file=xgif)
	with open(outygi, "w") as ygif:
		for key in peaks:
			print(key, file=ygif)
	end_time = pc()
	print('Used (secs): ', end_time - start_time)

if __name__ == "__main__":
	"""count reads to peak regions"""
	run()


