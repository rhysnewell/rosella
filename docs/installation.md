---
title: Installation
---

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
```
conda env create -n rosella -f rosella.yml
```

## Shell completion

Completion scripts for various shells e.g. BASH can be generated. For example, to install the bash completion script system-wide (this requires root privileges):

```
rosella shell-completion --output-file rosella --shell bash
mv rosella /etc/bash_completion.d/
```

It can also be installed into a user's home directory (root privileges not required):

```
rosella shell-completion --shell bash --output-file /dev/stdout >>~/.bash_completion
```

In both cases, to take effect, the terminal will likely need to be restarted. To test, type `rosella rec` and it should complete after pressing the TAB key.