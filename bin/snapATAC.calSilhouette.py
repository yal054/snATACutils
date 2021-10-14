#!/bin/python

import argparse

parser = argparse.ArgumentParser(description='calculate silhouette score for each clustering')
parser.add_argument('-i', '--meta', type=str, dest="meta", help='input meta table with umap')
parser.add_argument('-c', '--cluster', type=str, dest="cluster", help='partition from leiden')
parser.add_argument('-o', '--out', type=str, dest="output", help='prefix of output file')

args = parser.parse_args()

from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
from time import perf_counter as pc
import warnings

def run():
    warnings.simplefilter("ignore")
    """ Run and count to peak """
    start_time = pc()
    """ init input files """
    metaf = args.meta
    clsf = args.cluster
    outf = args.output
    
    # load and parse data
    meta = pd.read_csv(metaf, sep="\t", header=0)
    cluster_labels = np.loadtxt(clsf)
    X = meta[['umap-1', 'umap-2']].to_numpy()
    
    # the average value for all the cells.
    n_clusters = int(max(cluster_labels))
    silhouette_avg = silhouette_score(X, cluster_labels)
    print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg)

    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    # plot 
    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)
    
    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])
    y_lower = 10
    
    centers = np.empty((0,2))
    for i in range(n_clusters):
        i = i + 1
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]
    
        ith_cluster_silhouette_values.sort()
    
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
    
        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                              0, ith_cluster_silhouette_values,
                              facecolor=color, edgecolor=color, alpha=0.7)
    
        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
    
        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples
        
        # centers of each clusters
        center = X[cluster_labels == i]
        center = np.nanmean(center,axis=0)
        centers = np.vstack([centers, center])
    
    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")
    
    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")
    
    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    
    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(X[:, 0], X[:, 1], marker='.', s=30, lw=0, alpha=0.7,
                    c=colors, edgecolor='k')
    
    # Draw white circles at cluster centers
    ax2.scatter(centers[:, 0], centers[:, 1], marker='o',
                    c="white", alpha=1, s=200, edgecolor='k')
    
    for i, c in enumerate(centers):
        ax2.scatter(c[0], c[1], marker='$%d$' % i, alpha=1,
                        s=50, edgecolor='k')
    
    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("UMAP-1")
    ax2.set_ylabel("UMAP-2")
    
    plt.suptitle(("Silhouette analysis for KMeans clustering on sample data "
                      "with n_clusters = %d" % n_clusters),
                     fontsize=14, fontweight='bold')
    
    fileN = [outf, "silhouette", "png"]
    fileN = '.'.join(fileN)
    fig.savefig(fileN)

    outpart = ".".join([outf,"silhouette.txt"])
    with open(outpart, 'w') as filehandle:
        for listitem in sample_silhouette_values:
            listitem = listitem
            filehandle.write('%s\n' % listitem)
    end_time = pc()
    print('Used (secs): ', end_time - start_time)

if __name__ == "__main__":
    """calculate silhouette"""
    run()

