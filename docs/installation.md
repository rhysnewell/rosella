# Installation

Rosella is one binary with nothing to install beside it. It maps short reads with strobealign
0.19 and long reads with minimap2 2.31, and searches and aligns genes with tongs, a Rust HMMER
3.4, all built in. It calls the genes itself.

## Conda

```bash
conda create -n rosella -c bioconda rosella
conda activate rosella
rosella --version
```

## Cargo

With Rust 1.88 or later and a C compiler, from crates.io:

```bash
cargo install --locked rosella
```

Or straight from GitHub:

```bash
cargo install --locked --git https://github.com/rhysnewell/rosella
```

`--locked` builds the dependency versions rosella was tested with.

## macOS binary

Each [release](https://github.com/rhysnewell/rosella/releases) ships a binary for arm64 macOS.
Unpack it and put `rosella` on your `PATH`. On Linux, install with conda or cargo.

## The GPU

`--gpu` trains the learned space, embeds the contigs and runs the big neighbour searches on the
GPU. A run that cannot use one logs why on a line starting `GPU: none` and runs on the CPU.

On macOS it uses Metal. Metal ships with the system, so nothing else is needed.

On Linux it uses CUDA. Rosella loads the CUDA driver, cuBLAS and NVRTC from CUDA 12 when the run
starts. A build needs none of them, but a `--gpu` run needs all three on the library path. On a
cluster that is often a `module load CUDA` before the run.

## Checking the install

```bash
rosella --version
```

The version carries the commit it was built from, and a `-dirty` suffix when the tree it was
built in had uncommitted changes. Quote that string, not the release number, when you report a
result or a bug.
