# Installation

Rosella is one static binary plus the external tools it shells out to.

| Tool | When it is needed |
|---|---|
| `hmmalign`, from HMMER 3.4 or newer | Only for `--checkm` and `rosella score`, to align duplicated CheckM markers |
| `minimap2` and `samtools` | Only to map reads. BAMs and `--coverage-file` need neither |
| `coverm` 0.6.1 or newer | Only to map reads with a mapper other than minimap2 |

Rosella calls the genes itself and searches them with HMMER 3.4 built in, so there is no gene
caller and no `hmmsearch` to install.

## Conda

```bash
conda create -n rosella -c bioconda rosella
conda activate rosella
rosella --version
```

## From source

With `rust` and `cargo` on your system:

```bash
git clone https://github.com/rhysnewell/rosella
cd rosella
cargo install --path .
```

That builds the binary but not the tools beside it. Get those from the environment file in
the repository:

```bash
mamba env create -f rosella.yml -n rosella
mamba activate rosella
rosella --help
```

## With pixi

```bash
git clone https://github.com/rhysnewell/rosella
cd rosella
cargo install --path .
pixi shell
rosella --help
```

## Checking the install

```bash
rosella --version
```

The version carries the commit it was built from, and a `-dirty` suffix when the tree it was
built in had uncommitted changes. Quote that string, not the release number, when you report a
result or a bug.
