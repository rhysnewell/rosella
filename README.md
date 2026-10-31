[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/rosella/README.html)
![](https://anaconda.org/bioconda/rosella/badges/license.svg)
![](https://anaconda.org/bioconda/rosella/badges/version.svg)
![](https://anaconda.org/bioconda/rosella/badges/platforms.svg)
[![DOI](https://zenodo.org/badge/247065826.svg)](https://zenodo.org/doi/10.5281/zenodo.10140531)

![Rosella logo](images/rosella.png)

# Rosella

Rosella recovers genomes from a metagenome assembly. It learns a space from each contig's
composition and depth, and joins every contig to its nearest neighbours in that space. The
neighbour graph is cut with Leiden and label propagation across a ladder of resolutions.
Single copy markers pick which cut each bin takes. Bins that miss the bar go back into a pool to
be searched again.

All of it is one binary, written in Rust. Reads are mapped in process with strobealign 0.19 and
minimap2 2.31. BAMs are read in one pass for each contig's depth and the depth of every 100 bp
window along it. The learned space trains on those windows. Genes are called in process and searched
for the markers with tongs, a Rust HMMER 3.4. There is nothing else to install.

`--gpu` trains the learned space and runs the big neighbour searches on the GPU, with Metal on
macOS and CUDA on Linux. Without a GPU the same run goes on the CPU.

```bash
rosella recover -r assembly.fasta -b sample_*.bam -o rosella_bins/ -t 24 --gpu
```

## Install

### Conda

```bash
conda create -n rosella -c bioconda rosella
conda activate rosella
rosella --version
```

### Cargo

With Rust 1.88 or later and a C compiler, from crates.io or straight from GitHub:

```bash
cargo install --locked rosella
cargo install --locked --git https://github.com/rhysnewell/rosella
```

### macOS binary

Each [release](https://github.com/rhysnewell/rosella/releases) ships a binary for arm64 macOS.

## Documentation

The [documentation](https://rhysnewell.github.io/rosella) covers each input, the output and how
the binner works.

## License

Rosella is licensed under the GNU General Public License v3.0 only. See [LICENSE](LICENSE).

Copyright (c) 2020, Centre for Microbiome Research, QUT
