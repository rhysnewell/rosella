#!/usr/bin/env python3

import umap
from threadpoolctl import threadpool_limits
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyfaidx
import sys
from sklearn.cluster import OPTICS

#random.seed(a=345210)

def bin_contigs(bin_dict, assembly_file):
    assembly = pyfaidx.Faidx(assembly_file)
    for (bin, contigs) in bin_dict.items():
        with open('bin.'+str(bin)+'.fna', 'w') as f:
            for contig in contigs:
                seq = assembly.fetch(contig, 1, assembly.index[contig].rlen)
                fasta = ">" + seq.name + '\n'
                fasta += seq.seq + '\n'
                f.write(fasta)

if __name__=="__main__":
    try:
        coverages = sys.argv[1]
        fasta = sys.argv[2]
        min_dist = float(sys.argv[3])
        spread = float(sys.argv[4])
        n_neighbours = int(sys.argv[5])
        max_eps = float(sys.argv[6])

    except IndexError:
        print("Usage <Coverages> <Fasta> <min_dist> <spread> <n_neighbours> <epsilon>")
        sys.exit()

    # Contig coverage array
    coverage_array = []
    # Idx : Contig Name
    coverage_key = {}
    with open(coverages) as file:
        for (idx, line) in enumerate(file):
            if idx == 0:
                continue
            else:
                line = line.split()
                coverage_key[idx] = line[0]
                coverage_array.append([float(val) for val in line[3:]])

    coverage_array = np.array(coverage_array)
    reducer = umap.UMAP(metric='manhattan', min_dist=min_dist, spread=spread, n_neighbors=n_neighbours)
    with threadpool_limits(limits=20, user_api='blas'):
        coverage_array = coverage_array / np.linalg.norm(coverage_array)
        embedding = reducer.fit_transform(coverage_array)
        clustering = OPTICS(max_eps=max_eps).fit(embedding)
        plt.scatter(embedding[:, 0], embedding[:, 1], c=clustering.labels_)
        plt.gca().set_aspect('equal', 'datalim')
        plt.title('UMAP projection of contig clusters', fontsize=24)
        plt.savefig("umap_projection.png")

        contig_clusters = {}
        for (idx, label) in enumerate(clustering.labels_):
            idx += 1
            try:
                contig_clusters[label].append(coverage_key[idx])
            except KeyError:
                contig_clusters[label] = [coverage_key[idx]]

        bin_contigs(contig_clusters, fasta)
