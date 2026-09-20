# Usage

## Coverage you already have

The common case. You are running rosella beside other binners and already have a CoverM contig
table, made with `coverm contig -m metabat`:

```bash
rosella recover -r assembly.fasta -C coverage.tsv -o rosella_bins/ -t 24
```

## Let rosella map the reads

Pass reads instead and rosella hands them to CoverM. Long and short reads can go in together,
and so can as many samples as you have. More samples is better: coverage over one sample is a
single number per contig, and a lot of what the binner can separate comes from how that number
moves across samples.

```bash
rosella recover -r assembly.fasta \
  -1 sample_{1,2,3}.1.fq.gz -2 sample_{1,2,3}.2.fq.gz \
  --longreads nanopore.fq.gz \
  -o rosella_bins/ -t 24
```

Sorted BAM files work too, with `-b` for short read and `-l` for long read, and skip the mapping.

Rosella does not keep the BAMs it maps into. The coverage it computed lands in the output
directory as `coverage.tsv`, and a later run pointed at that directory picks it up on its own or
takes it back with `-C`. The composition table is only kept with `--write-kmer-table`, since it
runs to hundreds of megabytes on a large assembly and costs seconds to rebuild.

## Refining bins that already exist

`rosella refine` takes bins from any binner and re-partitions each one on its own. It needs the
assembly the bins were built from, the bins, and coverage:

```bash
rosella refine -r assembly.fasta -d metabat_bins/ -x fa -C coverage.tsv -o refined_bins/ -t 24
```

Pass a bin quality table with `--bin-quality` to aim the work: a CheckM1, CheckM2 or AMBER table
is read as it is, and only bins over `--split-contamination`, 15 per cent by default, are always
treated as split candidates.

`--output-directory` has to hold no `.fna` files already. A second run into a used directory is
refused rather than appended to. Bins the refiner did not change are written out beside the ones
it split, renumbered, so the output names do not carry over from the input.

## Scoring bins

`rosella score` runs the marker annotation over a set of bins and writes a completeness and
contamination table. No gold standard and no reference database:

```bash
rosella score -r assembly.fasta -d rosella_bins/ -x fna -o quality.tsv -t 24
```

## Reusing the marker annotation

Calling and searching the genes is most of the wall time of a run. `--marker-cache <dir>` keys
the annotation on the assembly and on every setting that changes it, so a second run over the
same assembly skips it:

```bash
rosella recover -r assembly.fasta -C coverage.tsv -o rosella_bins/ --marker-cache ~/.cache/rosella
```

The key includes the gene caller's own version, so upgrading rosella can correctly invalidate a
cached annotation and pay for it once.

## Reproducibility

Two runs of one assembly with the same inputs give the same bins. The partition ensemble and the
neighbour search take fixed seeds rather than `--seed`, so there is nothing to choose and no
reason to run three times. `--seed`, `--knn-seed` and `--partition-seed` exist to probe that, not
to be tuned.

## Flags worth knowing

| Flag | Default | What it does |
|---|---|---|
| `--min-contig-size` | 1500 | Contigs shorter than this take no part in binning |
| `--min-bin-size` | 200000 | Clusters smaller than this are not written as a bin |
| `-t, --threads` | 10 | Threads for rosella and for everything it calls |
| `--no-refine` | off | Skip the splitter, which cuts up the chimeric bins the first clustering leaves behind |
| `--recruit` | off | Let a bin short of the bars take single contigs back off its neighbours |
| `--assembly-graph` | none | A GFA whose links join the neighbour graph as extra edges |

Everything else is in `rosella recover --full-help`, grouped by what it touches. The flags under
"Reports" change no bin and only write a table about the run.
