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


def write_contig(contig, assembly, f):
    seq = assembly.fetch(contig, 1, assembly.index[contig].rlen)
    fasta = ">" + seq.name + '\n'
    fasta += seq.seq + '\n'
    f.write(fasta)


def bin_contigs(bin_dict, assembly_file, min_bin_size=200000):
    assembly = pyfaidx.Faidx(assembly_file)
    for (bin, contigs) in bin_dict.items():
        if bin != -1:
            # Calculate total bin size and check if it is larger than min_bin_size
            bin_length = sum([assembly.index[contig].rlen for contig in contigs])
            if bin_length >= min_bin_size:
                with open('bin.'+str(bin)+'.fna', 'w') as f:
                    for contig in contigs:
                        write_contig(contig, assembly, f)

        else:
            # Get final bin value
            bin_max = max(bin_dict.keys()) + 1
            # Rescue any large unbinned contigs and put them in their own cluster
            for contig in contigs:
                if assembly.index[contig].rlen >= min_bin_size:
                    with open('bin.' + str(bin_max) + '.fna', 'w') as f:
                        write_contig(contig, assembly, f)
                    bin_max += 1

if __name__=="__main__":
    try:
        coverages = sys.argv[1]
        fasta = sys.argv[2]
        min_dist = float(sys.argv[3])
        spread = float(sys.argv[4])
        n_neighbours = int(sys.argv[5])
        max_eps = float(sys.argv[6])
        min_length = int(sys.argv[7])
        min_cls_size = 200000

    except IndexError:
        print("Usage <Coverages> <Fasta> <min_dist> <spread> <n_neighbours> <epsilon> <min contig length>")
        sys.exit()

    # Contig coverage array
    coverage_array = []
    # Idx : Contig Name
    coverage_key = {}
    with open(coverages) as file:
        index = 1
        for (idx, line) in enumerate(file):
            if idx == 0:
                continue
            else:
                line = line.strip().split()
                if float(line[1]) >= min_length:
                    coverage_key[index] = line[0]
                    coverage_array.append([float(val) for val in line[3:]])
                    index += 1

    coverage_array = np.array(coverage_array)
    reducer = umap.UMAP(metric='manhattan', min_dist=min_dist, spread=spread, n_neighbors=n_neighbours)
    with threadpool_limits(limits=20, user_api='blas'):
        coverage_array = coverage_array / np.linalg.norm(coverage_array)
        embedding = reducer.fit_transform(coverage_array)
        clustering = OPTICS(min_samples=2).fit(embedding)
        plt.scatter(embedding[:, 0], embedding[:, 1], c=clustering.labels_, s=5, alpha=0.33)
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

