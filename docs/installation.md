# Installation

Rosella is one static binary with nothing to install beside it. It maps short reads with
strobealign 0.19 and long reads with minimap2 2.31, and searches and aligns genes with tongs, a Rust HMMER 3.4,
all built in. It calls the genes itself.

## Conda

```bash
conda create -n rosella -c bioconda rosella
conda activate rosella
rosella --version
```

## From source

With `rust`, `cargo` and a C compiler on your system:

```bash
git clone https://github.com/rhysnewell/rosella
cd rosella
cargo install --path .
```

## Checking the install

```bash
rosella --version
```

The version carries the commit it was built from, and a `-dirty` suffix when the tree it was
built in had uncommitted changes. Quote that string, not the release number, when you report a
result or a bug.
