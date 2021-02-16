#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='leiden on 3-cols knn sparse matrix')
parser.add_argument('-i', '--in', type=str, dest="input", help='input 3-col file: peak, cell, ct')
parser.add_argument('-r', '--resolution', type=float, dest="resolution", help='resolution from 0-1')
parser.add_argument('-o', '--out', type=str, dest="output", help='prefix of output file')

args = parser.parse_args()

import leidenalg as la
import igraph as ig
from scipy.io import mmread
from scipy import sparse
from time import perf_counter as pc

def run():
    """ Run and count to peak """
    start_time = pc()
    """ init input files """
    inf = args.input
    reso = args.resolution
    outf = args.output
    knn = mmread(inf)
    vcount = max(knn.shape)
    sources, targets = knn.nonzero()
    edgelist = list(zip(sources.tolist(), targets.tolist()))
    g = ig.Graph(vcount, edgelist)
    
    partition = la.find_partition(g, la.RBConfigurationVertexPartition, resolution_parameter = reso)
    #partition = la.find_partition(g, la.CPMVertexPartition, resolution_parameter = reso)
    #partition = la.find_partition(g, la.RBERVertexPartition, resolution_parameter = reso)
    part_membership = partition.membership

    outpart = ".".join([outf,"partition.txt"])
    with open(outpart, 'w') as filehandle:
        for listitem in part_membership:
            listitem = listitem + 1
            filehandle.write('%s\n' % listitem)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

if __name__ == "__main__":
    """leiden on knn"""
    run()


