#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='leiden on 3-cols knn sparse matrix')
parser.add_argument('-i', '--in', type=str, dest="input", help='input 3-col file: peak, cell, ct')
parser.add_argument('-r', '--resolution', type=float, dest="resolution", help='resolution from 0-1')
parser.add_argument('--u1', type=float, dest="u1", default=0.05, help='left interval cutoff of CDF')
parser.add_argument('--u2', type=float, dest="u2", default=0.95, help='right interval cutoff of CDF')
parser.add_argument('-n', type=int, dest="N", help='iteration of N time')
#parser.add_argument('-c', '--cpu', type=int, dest="cpu", default=1, help='# of cpu')
parser.add_argument('-o', '--out', type=str, dest="output", help='prefix of output file')

args = parser.parse_args()

#from multiprocessing import Pool
import leidenalg as la
import igraph as ig
from scipy.io import mmread
from scipy import sparse
from time import perf_counter as pc
from scipy.sparse import csr_matrix
from scipy.sparse import save_npz
import numpy as np
import itertools
from scipy.cluster.hierarchy import linkage, leaves_list, cophenet
from scipy.spatial.distance import squareform

import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig, imshow, set_cmap

plt.switch_backend('agg')
#try:
#    from matplotlib.pyplot import savefig, imshow, set_cmap
#except ImportError as exc:
#    warn("Matplotlib must be installed.")
import fastcluster as fc

def cal_connectivity(P):
    """ calculate connectivity matrix """
    #print("=== calculate connectivity matrix ===")
    #connectivity_mat = csr_matrix((len(P), len(P)))
    connectivity_mat = np.zeros((len(P), len(P)))
    classN = max(P)
    for cls in range(classN + 1):
        xidx = [i for i, x in enumerate(P) if x == cls]
        iterables = [ xidx, xidx ]
        for t in itertools.product(*iterables):
            connectivity_mat[t[0],t[1]] = 1
    """connectivity_mat = csr_matrix(connectivity_mat)"""
    return connectivity_mat

def plotC(prefix, C):
    """
    Plot reordered consensus matrix.

    :param C: Reordered consensus matrix.
    :type C: numpy.ndarray`
    :param rank: Factorization rank.
    :type rank: `int`
    """
    fig = plt.figure(figsize=(5,5), dpi=100);
    imshow(C)
    set_cmap("RdBu_r")
    fileN = [prefix, "C", "png"]
    fileN = '.'.join(fileN)
    fig.savefig(fileN)

def reorder(C):
    """
    Reorder consensus matrix.

    :param C: Consensus matrix.
    :type C: `numpy.ndarray`
    """
    Y = 1 - C
    Z = linkage(squareform(Y), method='average')
    ivl = leaves_list(Z)
    ivl = ivl[::-1]
    return C[:, ivl][ivl, :]


def plot_CDF(prefix, C, u1, u2, num_bins=100):
    counts, bin_edges = np.histogram(C, bins=num_bins, density=True)
    cdf = np.cumsum(counts)
    fig = plt.figure(dpi=100);
    plt.plot(bin_edges[1:], cdf)
    plt.xlabel('Consensus index value')
    plt.ylabel('CDF')
    plt.axvline(x=u1, color='grey', linestyle='--')
    plt.axvline(x=u2, color='grey', linestyle='--')
    fileN = [prefix, "cdf", "png"]
    fileN = '.'.join(fileN)
    fig.savefig(fileN)
    outBinEdges = ".".join([prefix, "cdf.txt"])
    with open(outBinEdges, 'w') as fo:
        fo.write('\t'.join(str(i) for i in cdf) + "\n")


def cumfreq(a, numbins=100, defaultreallimits=None):
    # docstring omitted
    h,l,b,e = histogram(a,numbins,defaultreallimits)
    cumhist = np.cumsum(h*1, axis=0)
    return cumhist,l,b,e


def cal_cophenetic(C):
    """ calculate cophenetic correlation coefficient """
    print("=== calculate cophenetic correlation coefficient ===")
    X = C
    """Z = linkage(X)"""
    Z = fc.linkage_vector(X)         # Clustering
    orign_dists = fc.pdist(X)  # Matrix of original distances between observations
    cophe_dists = cophenet(Z)  # Matrix of cophenetic distances between observations
    corr_coef = np.corrcoef(orign_dists, cophe_dists)[0,1]
    return corr_coef

def cal_dispersion(C):
    """ calculate dispersion coefficient """
    print("=== calculate dispersion coefficient ===")
    start_t = pc()
    n = C.shape[1]
    corr_disp = np.sum(4 * np.square(np.concatenate(C - 1/2)))/(np.square(n))
    end_t = pc()
    print('Used (secs): ', end_t - start_t)
    return corr_disp

def cal_PAC(C, u1, u2):
    """ calculate PAC (proportion of ambiguous clustering) """
    print("=== calculate PAC (proportion of ambiguous clustering) ===")
    start_t = pc()
    totalN = C.shape[0] * C.shape[0]
    u1_fraction = (C.ravel() < u1).sum() / totalN
    u2_fraction = (C.ravel() < u2).sum() / totalN
    PAC = u2_fraction - u1_fraction
    end_t = pc()
    print('Used (secs): ', end_t - start_t)
    return PAC

def cal_stab(x):
    """ calculate stability for every cell """
    s = np.sum(abs(x - 0.5)) / (0.5 * x.shape[0])
    return s

def run():
    """ Run and count to peak """
    start_time = pc()
    """ init input files """
    inf = args.input
    reso = args.resolution
    N = args.N
#    ncpu = args.cpu
    u1 = args.u1
    u2 = args.u2
    outf = args.output

    knn = mmread(inf)
    knn = knn.tocsr()
    dimN = knn.shape[0]
    if dimN > 5000:
        import random
        random.seed(2020)
        idxy = random.sample(range(dimN), 5000)
        knn = knn[idxy, :]
        knn = knn[:, idxy] # downsample (5,000 observations)
    vcount = max(knn.shape)
    sources, targets = knn.nonzero()
    edgelist = list(zip(sources.tolist(), targets.tolist()))
    g = ig.Graph(vcount, edgelist)

    # generate consensus matrix    
    dims = knn.shape[0]
    consensus = csr_matrix((dims, dims))
    print("=== calculate connectivity matrix ===")
    for seed in range(N):
        start_t =pc()
        partition = la.find_partition(g, la.RBConfigurationVertexPartition, resolution_parameter = reso, seed = seed)
        part_membership = partition.membership
        consensus += cal_connectivity(part_membership)
        end_t = pc()
        print('seed ', seed, ' used (secs): ', end_t - start_t)
    consensus /= N

    # save consensus matrix
    consensus_sp = sparse.csr_matrix(consensus)
    outfname = '.'.join([outf, "consensus", "npz"])
    save_npz(outfname, consensus_sp)

    # plotting
    order_consensus = reorder(consensus)
    plotC(outf, order_consensus)
    plot_CDF(outf, consensus, u1, u2)

    # cal measurement
#    o_cophcor = cal_cophenetic(consensus)
    o_disp = cal_dispersion(consensus)
    o_PAC = cal_PAC(consensus, u1, u2)
    print("=== calculate stability for every cell ===")
    o_stab = np.apply_along_axis(cal_stab, 0, consensus)

    # write stat
    out_list = [outf, reso, o_disp, o_PAC]
    outstat = ".".join([outf, "stat.txt"])
    with open(outstat, 'w') as fo:
        fo.write('\t'.join(str(i) for i in out_list) + '\n')

    outStab = ".".join([outf, "stab.txt"])
    with open(outStab, 'w') as fo:
        fo.write('\t'.join(str(i) for i in o_stab) + '\n')

    end_time = pc()
    print('Total used (secs): ', end_time - start_time)


if __name__ == "__main__":
    """consensus leiden on knn"""
    run()


