# Installation

Rosella is one static binary plus two external tools it shells out to.

| Tool | When it is needed |
|---|---|
| `hmmsearch`, from HMMER 3.4 or newer | Always. The single copy marker search |
| `coverm` 0.6.1 or newer | Only when you do not pass `--coverage-file` |

Rosella calls the genes itself, so there is no gene caller to install.

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

That builds the binary but not the two tools beside it. Get those from the environment file in
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
