# Output

Everything a run produces lands in `--output-directory`.

## Bins

| Name | What is in it |
|---|---|
| `rosella_bin_N.fna` | A bin. Numbered from zero |
| `rosella_bin_single_contig_N.fna` | One contig that clears `--min-bin-size` on its own |
| `rosella_bin_replicon_N.fna` | One contig of 10 kb or more that sat in a bin but carries no single copy marker, has the dense, short gene layout of a phage or plasmid, and is further from the bin's composition than any of its marker carrying contigs, or is a passenger with that layout. Published alone whatever its size |
| `rosella_bin_unbinned.fna` | Contigs that took part in binning and found no bin, contigs from `--attach-floor` to `--min-contig-size` that no bin took, and passengers: marker-free contigs of 10 kb or more that sat in a bin but lie beyond 98 per cent of the assembly's marker carrying contigs of like length on both composition and depth |
| `rosella_bin_small_unbinned.fna` | Contigs under `--attach-floor`, which never took part |

Every contig in the assembly is written exactly once, across those five. A run that cannot
account for all of them fails rather than writing a partial answer.

A `refine` run names its bins after `--bin-tag`, `refined_1` by default.

## Tables

| File | What it holds |
|---|---|
| `coverage.tsv` | The coverage table, only when rosella computed it. A table passed with `-C` is read where it is and not copied |
| `kmer_frequencies.k4.tsv` | The composition table, only with `--write-kmer-table` |
| `quality.tsv` | Completeness and contamination per bin, on rosella's markers and on CheckM1's |
| `timings.tsv` | Seconds and percent per stage |
| `stages.tsv` | Bin counts and sizes after each stage of the run |

A `coverage.tsv` or a composition table already sitting in the output directory is picked up by
a later run pointed at it, which is also why a run with different inputs needs a different one.

## The quality table

`quality.tsv` scores every bin on two marker panels. Both come from the one annotation the
binner already made, and the same lineage pick decides which set each panel reads.

| Column | What it holds |
|---|---|
| `set` | The lineage the bin's markers fit: `bac`, `ar`, `pat` (Patescibacteria) or `dpann` |
| `gtdb_completeness`, `gtdb_contamination` | The panel the binner decides with, built from GTDB's bac120 and ar53 markers and counted per marker |
| `checkm_completeness`, `checkm_contamination` | CheckM1's own marker sets, groups and arithmetic (Parks et al. 2015) |

The CheckM columns read Bacteria for `bac`, Archaea for `ar` and the 43 CPR markers of Brown et
al. (2015) for `pat`. CheckM ships no DPANN set, so `dpann` reads rosella's DPANN markers grouped
the way CheckM groups its own. Compare a bin with CheckM1 on these columns. On 723 bacterial and
Patescibacteria bins from two sludge assemblies they read within about a point of CheckM1 on
bacterial bins, and about a point under CheckM2 on contamination.

The GTDB panel holds genes that are rarely duplicated within a genome, so a second copy is strong
evidence of a second organism. That suits binning decisions. It also makes its contamination read
lower than CheckM's on the same bin.

Neither is a substitute for CheckM2 in a paper. Both are the binner's own view of its own work, and
a bin that looks clean to the markers can still carry foreign sequence the markers cannot see:
contamination in base pairs and contamination in duplicated marker genes are different numbers.

## Reports

The flags under "Reports" in `--full-help` each write one table about the run and change no bin,
so a run with one of them on is directly comparable to a run without it. `--knn-report` and
`--reach-report` stop before partitioning and write nothing else.
