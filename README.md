[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/rosella/README.html)
![](https://anaconda.org/bioconda/rosella/badges/license.svg)
![](https://anaconda.org/bioconda/rosella/badges/version.svg)
![](https://anaconda.org/bioconda/rosella/badges/platforms.svg)
[![DOI](https://zenodo.org/badge/247065826.svg)](https://zenodo.org/doi/10.5281/zenodo.10140531)

![Rosella logo](images/rosella.png)

# Rosella

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

Read mapping, gene calling and the marker search are built in, so nothing else is needed.

### macOS binary

Each [release](https://github.com/rhysnewell/rosella/releases) ships a binary for arm64 macOS.
Put it on your `PATH`:

```bash
install -m 755 rosella /usr/local/bin/
rosella --version
```

### GPU

On macOS `--gpu` runs on Metal. Metal ships with macOS. On Linux `--gpu` runs on CUDA. It needs
the NVIDIA driver plus cuBLAS and NVRTC from CUDA 12. Conda can install the two libraries:

```bash
conda create -n cuda -c conda-forge "cuda-version=12" libcublas cuda-nvrtc
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```

When no GPU can be opened, `--gpu` logs a warning and the run goes on the CPU.

## Examples

Map each sample's reads and keep their depth in one directory. `--interleaved` takes interleaved
reads, and `-b` takes BAMs already mapped:

```bash
rosella make -r assembly.fasta -o depth -t 16 \
    -1 s1_R1.fq.gz -2 s1_R2.fq.gz --longreads s2_ont.fq.gz
```

Bin from that depth:

```bash
rosella recover -r assembly.fasta --depth-directory depth -o bins -t 16
```

Refine bins from another binner:

```bash
rosella refine -r assembly.fasta --depth-directory depth -d other_bins -x fa -o refined -t 16
```

Pick the best bins from several binners:

```bash
rosella ensemble -r assembly.fasta --depth-directory depth -d bins_a bins_b -o ensemble -t 16
```

## Optional data

Rosella runs with what it ships and downloads nothing by default. Add `--download-data DIR` to any
run once to fetch CheckM's lineage sets, and any BUSCO sets asked for with `--auxiliary-markers`.
Later runs then report CheckM on the lineage set each bin's genome tree placement picks. The
directory is remembered inside the active conda or pixi environment, or under `~/.config/rosella`
outside one. `ROSELLA_DATA_DIR` points at a directory instead.

## License

Rosella is licensed under the GNU General Public License v3.0 only. See [LICENSE](LICENSE).

Copyright (c) 2020, Centre for Microbiome Research, QUT
