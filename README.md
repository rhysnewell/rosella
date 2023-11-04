![](https://travis-ci.com/rhysnewell/rosella.svg?branch=master)
![](https://anaconda.org/bioconda/rosella/badges/license.svg)
![](https://anaconda.org/bioconda/rosella/badges/version.svg)
![](https://anaconda.org/bioconda/rosella/badges/platforms.svg)
![Rosella logo](docs/_include/images/rosella.png)

# Rosella
Rosella is a metagenomic binning algorithm using UMAP and HDBSCAN. It is written in Rust with a python component that 
handles calls to UMAP and HDBSCAN. Rosella aims to be as user friendly as possible with multiple usage modes and installation
methods. 

Please note that Rosella is under active development with new commits often providing much improved results. If you would like
the most up to date version of Rosella please pull the code from `dev` branch. Hopefully releases will stabilise very soon.

## Quick Install
## Option 1: Conda

It's recommended that you create a new environment to ensure conda can correctly handle of the rosella's dependencies:

```bash
conda create -n rosella -c bioconda rosella
conda activate rosella
rosella --version
```

## Option 2: Install manually
After cloning the repo with `rust` and `cargo` installed on your system
```bash
cd rosella
cargo install --path .
```

Create the conda environment
```bash
mamba env create -f rosella.yml -n rosella
mamba activate rosella
rosella --help
```

## Documentation

Please refer to [documentation](https://rhysnewell.github.io/rosella) for installation and usage instructions.
