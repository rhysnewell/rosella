# Rosella
A metagenomic binner and variant clusterer using UMAP and HDBSCAN on compositonal datasets

# Installation

```
git clone https://github.com/rhysnewell/rosella.git
cd rosella
conda env create -n rosella -f rosella.yml
conda activate rosella
pip install --editable .
rosella bin --help
```

# Requirements

Initial requirements for binsnek can be downloaded using the `binsnek.yml`:
```
conda env create -n rosella -f rosella.yml
```

# Usage

To perform mag recovery:
```
rosella bin --assembly scaffolds.fasta --input coverm.cov --output output_dir/ --threads 24
```
