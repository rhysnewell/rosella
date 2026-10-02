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
| `quality.tsv` | Completeness and contamination per bin, on rosella's markers, and on CheckM1's with `--checkm`, with CheckM1's extended summary columns |
| `quality_sets.tsv` | Each bin read against every lineage set, not only the one it was read against in `quality.tsv` |
| `marker_duplicates.tsv` | Each marker a bin holds more than once, and the contigs that carry it |
| `abundance.tsv` | Each bin's mean depth and share of the binned depth, per sample |
| `timings.tsv` | Seconds and percent per stage |
| `stages.tsv` | Bin counts and sizes after each stage of the run |

A `coverage.tsv` or a composition table already sitting in the output directory is picked up by
a later run pointed at it, which is also why a run with different inputs needs a different one.

## The quality table

`quality.tsv` scores every bin on rosella's marker panel, and on CheckM1's as well when the run
is given `--checkm`. The same lineage pick decides which set each panel reads.

| Column | What it holds |
|---|---|
| `set` | The lineage the bin's markers fit: `bac`, `ar`, `pat` (Patescibacteria) or `dpann` |
| `gtdb_completeness`, `gtdb_contamination` | The panel the binner decides with, built from GTDB's bac120 and ar53 markers and counted per marker |
| `checkm_completeness`, `checkm_contamination` | CheckM1's own marker sets, groups and arithmetic (Parks et al. 2015). `NA` without `--checkm` |
| `checkm_strain_heterogeneity` | CheckM1's strain heterogeneity: the share of pairs of duplicated marker copies above 90 per cent amino acid identity, over the copies the CheckM columns count. High means the contamination is a second strain, low means a second organism. `NA` without `--checkm` |
| `gtdb_markers`, `gtdb_marker_groups` | Markers and marker groups in the set the bin was read against. The GTDB panel is flat, so the two agree |
| `gtdb_copies_0` to `gtdb_copies_5plus` | How many of those markers the bin holds no, one, two, three, four, or five or more copies of. A copy cut by a contig end counts toward one copy, never a second |
| `checkm_markers`, `checkm_marker_groups`, `checkm_copies_0` to `checkm_copies_5plus` | The same on CheckM1's set. `NA` without `--checkm` |
| `n50`, `longest_contig`, `mean_contig_length` | From the bin's contig lengths |
| `ambiguous_bases`, `gc`, `gc_std` | N bases, GC per cent, and the spread of GC about the bin's GC over contigs longer than 1 kbp, as CheckM1 reports them |
| `genes`, `coding_density`, `mean_gene_length` | Predicted genes, per cent of bases coding, and mean gene length in bp. `recover` calls no genes on the short contigs it attaches below the annotated floor, so these read over the contigs it did call. A coding density well under 75 per cent suggests eukaryotic sequence |

The CheckM columns read Bacteria for `bac`, Archaea for `ar` and the 43 CPR markers of Brown et
al. (2015) for `pat`. CheckM ships no DPANN set, so `dpann` reads rosella's DPANN markers grouped
the way CheckM groups its own. CheckM's models are searched once the bins are known, over the
contigs of genome bins alone, so they cost a run nothing without the flag. Replicons read `NA`.
Strain heterogeneity aligns each duplicated copy to its marker with `hmmalign`. A cached run holds
no proteins, so it calls genes again on the few contigs that carry a doubled CheckM marker. A
`dpann` bin's copies are GTDB hits, whose proteins are translated again from the gene's place.
`rosella score` always fills the CheckM columns, and reads the same as `recover --checkm` on the
same bins. Compare a bin with CheckM1 on these columns. On 723 bacterial and Patescibacteria bins
from two sludge assemblies they read within about a point of CheckM1 on bacterial bins, and about
a point under CheckM2 on contamination.

The GTDB panel holds genes that are rarely duplicated within a genome, so a second copy is strong
evidence of a second organism. That suits binning decisions. It also makes its contamination read
lower than CheckM's on the same bin.

Neither is a substitute for CheckM2 in a paper. Both are the binner's own view of its own work, and
a bin that looks clean to the markers can still carry foreign sequence the markers cannot see:
contamination in base pairs and contamination in duplicated marker genes are different numbers.

## The other quality tables

`rosella score` writes these two beside its `-o` table, and `recover` writes them on every run.
Each holds its CheckM rows only when the CheckM columns of `quality.tsv` are filled.

`quality_sets.tsv` reads every bin against every lineage set, one row per bin, set and panel, with
`chosen` set to 1 on the set `quality.tsv` reports. It shows how close the pick was, for example a
Patescibacteria bin read against `bac`.

`marker_duplicates.tsv` lists each marker of the chosen set that a bin holds more than once, one
row per contig carrying it, with its copies on that contig. A GTDB row gives each copy's place on
the contig as `begin-end:strand`. CheckM copies carry no place, so theirs read `NA`. A contig with
`copies` over 1 holds the marker twice on itself.

## Abundance

`abundance.tsv` has one row per bin of `quality.tsv`. For each sample of the coverage table it gives
the bin's mean depth weighted by contig length, and its share, in per cent, of depth times length
over every binned contig. `recover` only: `score` reads no coverage.

## Reports

The flags under "Reports" in `--full-help` each write one table about the run and change no bin,
so a run with one of them on is directly comparable to a run without it. `--knn-report` and
`--reach-report` stop before partitioning and write nothing else.

`--marker-report` writes one row per marker hit, with the bin its contig was written to first.
