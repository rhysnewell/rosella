[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/rosella/README.html)
![](https://anaconda.org/bioconda/rosella/badges/license.svg)
![](https://anaconda.org/bioconda/rosella/badges/version.svg)
![](https://anaconda.org/bioconda/rosella/badges/platforms.svg)
[![DOI](https://zenodo.org/badge/247065826.svg)](https://zenodo.org/doi/10.5281/zenodo.10140531)

![Rosella logo](images/rosella.png)

# Rosella
Rosella recovers genomes from a metagenome assembly by partitioning a k-nearest-neighbour graph over contig composition and coverage. Contigs are joined into a fuzzy simplicial set, the graph is cut with Leiden and label propagation across a ladder of resolutions, single-copy markers arbitrate which cut each bin takes, and the bins that miss the bar go back into a pool to be searched again. It is written entirely in Rust, with no
Python component and no external binning dependency.

The single-copy markers are searched and aligned in process, with tongs, a Rust HMMER 3.4, built into the
binary. Reads are mapped in process too, with strobealign 0.19 and minimap2 2.31 built in, so there is nothing else
to install. Rosella reads BAMs itself, in one pass for the coverage table and the window depth its learned
space uses. The genes
are called in process, so there is no gene caller to install.

```bash
rosella recover -r assembly.fasta -C coverage.tsv -o rosella_bins/ -t 24
```

Rosella is under active development and its results move between commits.

## Quick Install
## Option 1: Conda

It's recommended that you create a new environment to ensure conda can correctly handle of the rosella's dependencies:

```bash
conda create -n rosella -c bioconda rosella
conda activate rosella
rosella --version
```

## Option 2: Install manually
With `rust`, `cargo` and a C compiler installed on your system
```bash
git clone https://github.com/rhysnewell/rosella
cd rosella
cargo install --path .
```

## Documentation

Please refer to [documentation](https://rhysnewell.github.io/rosella) for installation and usage instructions.


## License

Rosella is licensed under the GNU General Public License v3.0 only. See [LICENSE](LICENSE).

Copyright (c) 2020, Centre for Microbiome Research, QUT
