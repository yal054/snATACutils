#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='Run statistic for NMF using sklearn.')
parser.add_argument('-m', '--matrix', type=str, dest="matrix", help='matrix in npz format')
parser.add_argument('-x', '--xgi', type=str, dest="xgi", help='input xgi index')
parser.add_argument('-y', '--ygi', type=str, dest="ygi", help='input ygi index')
parser.add_argument('--basis', type=str, dest="basis", help='basis matrix W')
parser.add_argument('--coef', type=str, dest="coef", help='coefficient matrix H')
#parser.add_argument('--conn', type=str, dest="conn", help='connectivity matrix C')
parser.add_argument('-c', '--contribute', type=float, dest="contribute", default=0.3, help='an float cutoff for contribute of bins / cells')
parser.add_argument('-o', '--outPrefix', type=str, dest="outPrefix", help='output prefix')

args = parser.parse_args()

from os.path import dirname, abspath, join
from warnings import warn

import numpy as np
from numpy import array, ravel
import scipy as sp
from scipy import io
from scipy.cluster.hierarchy import linkage, leaves_list, cophenet
from scipy.spatial.distance import squareform
from scipy.stats import describe
from scipy.stats.mstats import mquantiles
from scipy.sparse import save_npz, load_npz

from sklearn.decomposition import NMF
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_samples, silhouette_score, pairwise_distances

import fastcluster as fc

from time import perf_counter as pc

import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.switch_backend('agg')
try:
    from matplotlib.pyplot import savefig, imshow, set_cmap
except ImportError as exc:
    warn("Matplotlib must be installed.")

import itertools
import math
import copy
np.set_printoptions(suppress=True)

def run():
    """ Run standard NMF on rank """
    start_time = pc()
    """ init input files """
    V = args.matrix
    xgiF = args.xgi
    ygiF = args.ygi
    basisF = args.basis
    coefF = args.coef
    """connF = args.conn"""
    contribute = args.contribute
    outPrefix = args.outPrefix
    """perplexity = args.perplexity"""
    inF = read_files(V, xgiF, ygiF, basisF, coefF)
    V = inF['V']
    xgi = inF['xgi']
    ygi = inF['ygi']
    W = inF['basis']
    H = inF['coef']
    rank = inF['rank']
    """C = inF['conn']"""
    print("normalize W & H")
    normW = norm_W(W)
    np.savetxt('.'.join([outPrefix, "normW"]), normW, fmt= "%g", delimiter="\t")
    normH = norm_H(H)
    np.savetxt('.'.join([outPrefix, "normH"]), normH, fmt= "%g", delimiter="\t")
    print("calculate stat...")
    o_cell_class = def_cell_class(normH, contribute)
    o_cell_sparseness = cal_cell_sparseness(normH)
    o_entropy = cal_entropy(normH)    
    o_stat_H = stat_H(o_cell_class, o_cell_sparseness, o_entropy, ygi, outPrefix)    
    o_region_class = def_region_class(normW)
    o_featureScore_kim = cal_featureScore_kim(W)
    o_stat_W = stat_W(o_region_class, o_featureScore_kim, xgi, outPrefix)
    
    print("gerenate files based on stat...")
    gerenate_files(o_stat_H, o_stat_W, rank, V, outPrefix)
    """
    print("draw silhouette & tsne plot")
    o_silhouette = cal_silhouette(normH,o_stat_H)
    X_dist = cal_pairwise_pearson(normH)
    X_transformed = cal_tSNE(X_dist, perplexity)
    plot_silhouette_tsne(o_silhouette, X_transformed, o_stat_H, rank, outPrefix)
    """
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

def read_files(npz, xgi, ygi, basis, coef):
    """
    Read snATAC data in npz format. The matrix's shape is ### (bins) x ### (samples).
    Read in xgi and ygi
    Read in basis (W) and coefficient (H) matrix
    """
    V = load_npz(npz)
    V = V.tocsr()
    xfname = xgi
    with open(xfname) as xf:
        xgi_content = xf.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
    xgi_content = [x.strip() for x in xgi_content]    
    yfname = ygi
    with open(yfname) as yf:
        ygi_content = yf.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
    ygi_content = [y.strip() for y in ygi_content]    
    basis = np.loadtxt(basis)
    coef = np.loadtxt(coef)
    rank = coef.shape[0]
    """conn = np.loadtxt(conn)"""
    return {'V':V, 'xgi':xgi_content, 'ygi':ygi_content, 'basis':basis, 'coef':coef, 'rank':rank}

def norm_H(H):
    normH = H/H.sum(axis=0, keepdims=True)
    return normH

def norm_W(W):
    normW = W/W.sum(axis=1, keepdims=True)
    return normW

def def_cell_class(normH, contribute):
    """ define class from H with contributes cutoff"""
    print("=== extract feature from H with contribute cutoff ===")
    index = list(range(normH.shape[1]))
    colmax = np.amax(normH, axis=0)
    colsum = np.nansum(normH, axis=0)
    contributes = colmax / colsum
    class0 = normH.argmax(axis=0)
    unclass = np.where(contributes < contribute)
    unclass = unclass[0]
    unclassN = len(unclass)
    unclassPct = unclassN/normH.shape[1]
    class1 = copy.deepcopy(class0)
    class1[unclass] = 100
    stat = np.column_stack((index,class0,class1,contributes))
    stat[:,[0,1,2]] = stat[:,[0,1,2]].astype(int)
    return {'class0':class0, 'class1':class1, 'contributes':contributes, 'unclass':unclass, 'pct':unclassPct, 'stat':stat}

def cal_cell_sparseness(normH):
    print("=== calculate sparseness for each cell ===")
    sparseness = np.sqrt(np.nansum(normH**2, axis=0))
    rank = normH.shape[0]
    a = 1 / np.sqrt(rank)
    norm_sparse = (sparseness - a) / (1 - a) 
    return norm_sparse

def cal_entropy(normH):
    k = normH.shape[0]
    n = normH.shape[1]
    e_list = []
    for i in range(n):
        h_ij_list = []
        for j in range(k):
            h_ij = normH[j,i]
            if h_ij != 0:
                tmp = h_ij * math.log(h_ij, 2)
            else:
                tmp = 0
            h_ij_list.append(tmp)
        h_ij_sum = -np.nansum(h_ij_list)
        e_list.append(h_ij_sum)
    normInfoGain = 1 - np.nansum(e_list) / (n * math.log(k,2))
    return {'normInfoGain':normInfoGain, 'e_list':e_list}

def stat_H(o_cell_class, o_cell_sparseness, o_entropy, ygi, prefix):
    stat = o_cell_class['stat']
    norm_sparse = o_cell_sparseness
    entropy = o_entropy['e_list']
    statH = np.column_stack((stat,norm_sparse,entropy))
    statH[statH <= 0 ] = 0
    statH[:,[0,1,2]] = statH[:,[0,1,2]].astype(int)
    statH = np.column_stack((ygi,statH))
    head = "\t".join(["ygi", "index", "class0", "class1", "contributes", "sparseness", "entropy"])
    statH_fname = prefix + "." + "statH"
    np.savetxt(statH_fname, statH, fmt= "%s", delimiter="\t", header=head)
    return statH

def gerenate_files(o_stat_H, o_stat_W, rank, V, prefix):
    for r in range(rank):
        print(r)
        n = r + 1
        idy = np.concatenate(np.where(o_stat_H[:,2].astype(float).astype(int) == r))
        selt_ygi = o_stat_H[idy,0]
        selt_V = V[:,idy]
        ygi_fname = prefix + "." + "metacell_" + str(n) + "." + "ygi"
        V_fname = prefix + "." + "metacell_" + str(n) + "." + "npz"
        np.savetxt(ygi_fname, selt_ygi, fmt= "%s", delimiter="\t")
        save_npz(V_fname, selt_V)
        idx = np.concatenate(np.where((o_stat_W[:,2].astype(float).astype(int) == r)
            & (o_stat_W[:,5].astype(float).astype(int) == 1)
            & (o_stat_W[:,6].astype(float).astype(int) == 1)))
        selt_feature = o_stat_W[idx,0]
        feature_fname = prefix + "." + "metacell_" + str(n) + "." + "xgi"
        np.savetxt(feature_fname, selt_feature, fmt= "%s", delimiter="\t")


def cal_rss_mse(W, H, V):
    """ Residual Sum of Squares (RSS) & Mean Square Error (MSE)"""
    print("=== calculate Residual Sum of Squares (RSS) & Mean Square Error (MSE) ===")
    residualSquare = np.square(W.dot(H) - V)
    rss = np.nansum(residualSquare)
    mse = np.nanmean(residualSquare)
    out = [rss, mse]
    return out

def cal_evar(rss, V):
    print("=== calculate evar ===")
    evar = 1 - ( rss / np.nansum(V.data**2))
    return evar

def cal_featureScore_kim(W):
    """ extract feature from W """
    print("=== calculate featureScore from W ===")
    k = W.shape[1]
    m = W.shape[0]
    fs_list = []
    for i in range(m):
        rowsum = np.nansum(W[i,])
        p_iq_x_list = []
        for q in range(k):
            p_iq = W[i,q] / rowsum
            if p_iq != 0:
                tmp = p_iq * math.log(p_iq,2)
            else:
                tmp = 0
            p_iq_x_list.append(tmp)
        fs = 1 + 1/math.log(k,2) * np.nansum(p_iq_x_list)
        fs_list.append(fs)
    med = np.nanmedian(fs_list)
    mad = np.nanmedian(np.absolute(fs_list - med))
    fs_cutoff = med + 2 * mad
    selt_fs_idx = np.concatenate(np.where(fs_list >= fs_cutoff))
    return {'fs':fs_list, 'med': med, 'mad':mad, 'fs_cutoff': fs_cutoff, 'selt_fs_idx': selt_fs_idx}


def def_region_class(normW):
    """ define class from W with median contribute as cutoff"""
    print("=== extract feature from W with median contribute as cutoff ===")
    rowmax = np.amax(normW, axis=1)
    rowsum = np.nansum(normW, axis=1)
    contributes = rowmax / rowsum
    class0 = normW.argmax(axis=1)
    med_contribute = np.nanmedian(normW)
    selt_med_idx = np.concatenate(np.where(contributes > med_contribute))
    return {'class0':class0, 'contributes':contributes, 'med_contribute': med_contribute, 'selt_med_idx':selt_med_idx}

def stat_W(o_region_class, o_featureScore_kim, xgi, prefix):
    m = o_region_class['class0'].shape[0]
    index = list(range(m))
    class0 = o_region_class['class0']
    contributes = o_region_class['contributes']
    fs = o_featureScore_kim['fs']
    selt_fs_list = np.array([0] * m)
    selt_fs_idx = o_featureScore_kim['selt_fs_idx']
    selt_fs_list[selt_fs_idx] = 1
    selt_med_list = np.array([0] * m)
    selt_med_idx = o_region_class['selt_med_idx']
    selt_med_list[selt_med_idx] = 1
    statW = np.column_stack((index, class0, contributes, fs, selt_fs_list, selt_med_list))
    head = "\t".join(["chr", "start", "end", "index", "class0", "contributes", "fs", "selt_fs_list", "selt_med_list"])
    statW_fname = prefix + "." + "statW"
    statW[:,[0,1,4,5]] = statW[:,[0,1,4,5]].astype(int)
    statW = np.column_stack((xgi, statW))
    np.savetxt(statW_fname, statW, fmt= "%s", delimiter="\t", header=head)
    return statW
    
def cal_silhouette(normH,o_stat_H):
    X = normH.T
    n_clusters = normH.shape[0]
    cluster_labels = o_stat_H[:,2].astype(float).astype(int)
    """
    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed clusters
    """
    silhouette_avg = silhouette_score(X, cluster_labels)
    print("For n_clusters =", n_clusters, "The average silhouette_score is :", silhouette_avg)
    """# Compute the silhouette scores for each sample"""
    sample_silhouette_values = silhouette_samples(X, cluster_labels)
    return {'silhouette':silhouette_avg, 'silhouette_values':sample_silhouette_values}

def cal_nmds(X):
    embedding = MDS(n_components=2, dissimilarity='precomputed', n_jobs=2, verbose=2)
    X_transformed = embedding.fit_transform(X)
    return X_transformed

def cal_pairwise_dist(normH):
    X = normH.T
    X_pairwise_dist = pairwise_distances(X)
    return X_pairwise_dist

def cal_pairwise_pearson(normH):
    X = normH.T
    X_corrcoef = np.corrcoef(X)
    X_corrcoef_dist = np.sqrt(2*(1-X_corrcoef))
    return X_corrcoef_dist

def cal_tSNE(X_dist, p):
    X_transformed = TSNE(n_components=2, perplexity=p, random_state=1, verbose=2, metric="precomputed").fit_transform(X_dist)
    return X_transformed

def plot_silhouette_tsne(o_silhouette, X_transformed, o_stat_H, rank, prefix):
    n_clusters = rank
    silhouette_avg = o_silhouette['silhouette']
    sample_silhouette_values = o_silhouette['silhouette_values']
    cluster_labels = o_stat_H[:,2].astype(float).astype(int)

    """# Create a subplot with 1 row and 2 columns"""
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)
    fig.set_dpi(300)
    """
    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    """
    ax1.set_xlim([-0.1, 1])
    """
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    """
    ax1.set_ylim([0, len(X_transformed) + (n_clusters + 1) * 10])

    y_lower = 10
    for i in range(n_clusters):
        """# Aggregate the silhouette scores for samples belonging to cluster i, and sort them """
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
            0, ith_cluster_silhouette_values,
            facecolor=color, edgecolor=color, alpha=0.7)
        """# Label the silhouette plots with their cluster numbers at the middle"""
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        """# Compute the new y_lower for next plot"""
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel(" ".join(["The silhouette coefficient values (mean:",str(silhouette_avg),")"]))
    ax1.set_ylabel("Cluster label")

    """# The vertical line for average silhouette score of all the values """
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    """# 2nd Plot showing the actual clusters formed"""
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(X_transformed[:, 0], X_transformed[:, 1], marker='.', s=30, lw=0, alpha=0.7, c=colors, edgecolor='k')

    """ lable non-classified cells"""

    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("Feature space for the 1st tsne")
    ax2.set_ylabel("Feature space for the 2nd tsne")
    plt.suptitle(("Silhouette analysis for tSNE clustering on coefficient matrix H "
        "with n_clusters = %d" % n_clusters),fontsize=14, fontweight='bold')
    fig.savefig('.'.join([prefix, "silhouette_tsne", "png"]))

def cal_cophenetic(C):
    """ calculate cophenetic correlation coefficient """
    print("=== calculate cophenetic correlation coefficient ===")
    X = C  # Original data (1000 observations)
    """Z = linkage(X)"""
    Z = fc.linkage_vector(X)   # Clustering
    orign_dists = fc.pdist(X)  # Matrix of original distances between observations
    cophe_dists = cophenet(Z)  # Matrix of cophenetic distances between observations
    corr_coef = np.corrcoef(orign_dists, cophe_dists)[0,1]
    return corr_coef

def cal_dispersion(C):
    """ calculate dispersion coefficient """
    print("=== calculate dispersion coefficient ===")
    n = C.shape[1]
    corr_disp = np.nansum(4 * np.square(np.concatenate(C - 1/2)))/(np.square(n))
    return corr_disp    
        
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


if __name__ == "__main__":
    """Run statistics for NMF and save H & W."""
    run()
