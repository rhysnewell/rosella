---
title: Installation
---

Installation
========

## Option 1: Conda
This is the easiest method:

```
conda install -c bioconda rosella
rosella --version
```

It's recommended that you create a new environment to ensure conda can correctly handle of the rosella's dependencies:

```
conda create -n rosella -c bioconda rosella
conda activate rosella
rosella --version
```

## Option 2: Install manually
You may need to manually set the paths for `C_INCLUDE_PATH`, `LIBRARY_PATH`, `LIBCLANG_PATH`, 
`LD_LIBRARY_PATH`, and `OPENSSL_DIR` to their corresponding
paths in your conda environment if they can't properly be found on your system.
```

git clone --recursive https://github.com/rhysnewell/rosella.git \ 
cd rosella \
conda env create -n rosella -f rosella.yml \ 
conda activate rosella \ 
bash install.sh # or e.g. `cargo run -- bin`
```

Depending on your local network configuration, you may have problems obtaining rosella via git.
If you see something like this you may be behind a proxy that blocks access to standard git:// port (9418).

```
$ git clone --recursive git://github.com/rhysnewell/rosella.git
Cloning into 'rosella'...
fatal: Unable to look up github.com (port 9418) (Name or service not known)
```

Luckily, thanks to this handy tip from the developer of [Freebayes](https://github.com/ekg/freebayes) we can work around it.
If you have access to https:// on port 443, then you can use this 'magic' command as a workaround to enable download of the submodules:

```
git config --global url.https://github.com/.insteadOf git://github.com/
```

## Requirements

Initial requirements for rosella can be downloaded using the `rosella.yml`:
```
conda env create -n rosella -f rosella.yml
```