# Rosella

Rosella recovers metagenome assembled genomes from an assembly, using contig composition and
contig coverage and nothing else. It is written entirely in Rust. There is no Python component
and no other binner underneath it.

```bash
rosella recover -r assembly.fasta -C coverage.tsv -o rosella_bins/ -t 24
```

That is the whole thing for the common case. The pages here cover the other input shapes, what
lands in the output directory, and what the binner is doing in between.

## What it needs

An assembly, and coverage for it. Coverage can be BAMs, reads for rosella to map with
its built-in strobealign or minimap2, or a table you already have. BAMs and reads also give window depth, which a table cannot. Single copy markers are found in process: rosella calls the genes
itself and searches them with tongs, a Rust HMMER 3.4, built into the binary.

## The five subcommands

| Command | What it does |
|---|---|
| `rosella recover` | Bin an assembly from scratch |
| `rosella refine` | Refine bins that already exist, from rosella or another binner, through recover's refine cycle |
| `rosella ensemble` | Bin contigs as recover does, with other binners' bins offered beside its own partition |
| `rosella score` | Score a set of bins against the single copy markers, with no gold standard |
| `rosella make` | Write a depth directory from BAMs or reads, for `recover` and `refine` to read |

Every flag is documented in the binary. `rosella recover -h` lists the flags in everyday use and
`rosella recover --full-help` lists all of them, the same for `refine`, `ensemble`, `score` and
`make`. Where this guide and the binary disagree, the binary is right.

## Status

Rosella is under active development and its results move between commits. Pin a version if you
are comparing runs.
