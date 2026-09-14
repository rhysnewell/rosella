Installation
========

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

## Requirements

Initial requirements for rosella can be downloaded using the `rosella.yml`:
```bash
conda env create -n rosella -f rosella.yml
```
