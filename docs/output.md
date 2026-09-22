# Output

Everything a run produces lands in `--output-directory`.

## Bins

| Name | What is in it |
|---|---|
| `rosella_bin_N.fna` | A bin. Numbered from zero |
| `rosella_bin_single_contig_N.fna` | One contig that clears `--min-bin-size` on its own |
| `rosella_bin_replicon_N.fna` | One contig of 10 kb or more that sat in a bin but carries no single copy marker, has the dense, short gene layout of a phage or plasmid, and is further from the bin's composition than any of its marker carrying contigs. Published alone whatever its size |
| `rosella_bin_unbinned.fna` | Contigs that took part in binning and found no bin |
| `rosella_bin_small_unbinned.fna` | Contigs under `--min-contig-size`, which never took part |

Every contig in the assembly is written exactly once, across those five. A run that cannot
account for all of them fails rather than writing a partial answer.

A `refine` run names its bins after `--bin-tag`, `refined_1` by default.

## Tables

| File | What it holds |
|---|---|
| `coverage.tsv` | The coverage table, only when rosella computed it. A table passed with `-C` is read where it is and not copied |
| `kmer_frequencies.k4.tsv` | The composition table, only with `--write-kmer-table` |
| `quality.tsv` | Completeness and contamination per bin, from the single copy markers |
| `timings.tsv` | Seconds and percent per stage |
| `stages.tsv` | Bin counts and sizes after each stage of the run |

A `coverage.tsv` or a composition table already sitting in the output directory is picked up by
a later run pointed at it, which is also why a run with different inputs needs a different one.

## The quality table

`quality.tsv` is written from the same marker annotation the binner used to make its decisions,
in a CheckM1 shaped layout. It costs nothing extra, because the annotation already exists by the
time the bins are written.

It is not a substitute for CheckM2 in a paper. It is the binner's own view of its own work, and
a bin that looks clean to the markers can still carry foreign sequence the markers cannot see:
contamination in base pairs and contamination in duplicated marker genes are different numbers.

## Reports

The flags under "Reports" in `--full-help` each write one table about the run and change no bin,
so a run with one of them on is directly comparable to a run without it. `--knn-report` and
`--reach-report` stop before partitioning and write nothing else.
