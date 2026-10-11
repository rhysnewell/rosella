# Usage

## Coverage you already have

The common case. You are running rosella beside other binners and already have a CoverM contig
table, made with `coverm contig -m metabat`:

```bash
rosella recover -r assembly.fasta -C coverage.tsv -o rosella_bins/ -t 24
```

## Let rosella map the reads

Pass reads instead and rosella maps them with CoverM, with `minimap2 -x sr` for short reads and
`map-ont` for long reads unless `--mapper` or `--longread-mapper` says otherwise. Long and short
reads can go in together, and so can as many samples as you have. More samples is better: coverage over one sample is a
single number per contig, and a lot of what the binner can separate comes from how that number
moves across samples.

```bash
rosella recover -r assembly.fasta \
  -1 sample_{1,2,3}.1.fq.gz -2 sample_{1,2,3}.2.fq.gz \
  --longreads nanopore.fq.gz \
  -o rosella_bins/ -t 24
```

Each sample's sorted BAM is kept in `bams/` under the output directory, so a later run can take
them back with `-b`.

## BAMs you already have

Coordinate sorted BAMs skip the mapping, short and long reads alike, such as the ones `coverm
make` or aviary keeps:

```bash
rosella recover -r assembly.fasta -b binning_bams/*.bam -o rosella_bins/ -t 24
```

Rosella reads each BAM once. Every primary mapped alignment counts over its reference span, with
no identity filter and no end trimming. That one pass gives the coverage table and the depth of
every 100 bp window, and the learned space reads each view of a contig from those windows. A
coverage table alone carries no windows, so the space learns instead from the near complete bins
of a first pass in k-mer space, which recovers fewer genomes. Give BAMs or reads when you have
them.

The depth read from each BAM is kept in `depth/` under the output directory, or wherever
`--depth-directory` points:

| File | What it holds |
|---|---|
| `contigs.tsv` | Each contig and its length, in the BAMs' header order |
| `<sample>.win.f16` | The mean depth of every 100 bp window, as little-endian half floats, contig by contig in that order. A contig's last window covers what is left of it |
| `<sample>.depth.tsv` | That sample's whole-contig mean and variance, in CoverM's contig table layout |

A later run given the same BAMs and the same `--depth-directory` reads the kept depth and skips
the BAM, for every sample whose files are newer than its BAM. Other tools can read the windows
too, so binners run beside rosella need not walk the BAMs again.

## Depth without keeping the BAMs

`rosella make` writes the depth directory and nothing else, so the BAMs can be deleted as soon as
it finishes. It can run once per sample, each run adding to the directory, and `recover` reads
every sample in it in name order:

```bash
rosella make -b binning_bams/*.bam -o depth/ -t 24
rm binning_bams/*.bam
rosella recover -r assembly.fasta --depth-directory depth/ -o rosella_bins/ -t 24
```

Reads work too, with `-r` for the assembly to map them to.

The composition table is only kept with `--write-kmer-table`, since it runs to hundreds of
megabytes on a large assembly and costs seconds to rebuild.

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

`rosella score` calls genes on every contig in the bins it is given, searches them for the single
copy markers and writes a completeness and contamination table. No assembly, no gold standard
and no reference database:

```bash
rosella score -d rosella_bins/ -x fna -o quality.tsv -t 24
```

## Reusing the marker annotation

Calling and searching the genes is most of the wall time of a run. `--marker-cache <dir>` keys
the annotation on the assembly, its size and when it was last written, and on every setting that
changes it, so a second run over the same assembly skips it:

```bash
rosella recover -r assembly.fasta -C coverage.tsv -o rosella_bins/ --marker-cache ~/.cache/rosella
```

The key includes the gene caller's own version, so upgrading rosella can correctly invalidate a
cached annotation and pay for it once. The cache also keeps the CheckM search of every contig a
`--checkm` run has binned, so a later run over the same bins searches nothing again.

## Reproducibility

Two runs of one assembly with the same inputs give the same bins. The partition ensemble and the
neighbour search take fixed seeds rather than `--seed`, so there is nothing to choose and no
reason to run three times. `--seed`, `--knn-seed` and `--partition-seed` exist to probe that, not
to be tuned.

## Flags worth knowing

| Flag | Default | What it does |
|---|---|---|
| `--min-contig-size` | auto | Contigs at or above this are partitioned. auto reads 1500. A shorter cutoff takes nothing attach does not, and costs wall |
| `--attach-floor` | auto | Contigs from this length up to `--min-contig-size` join a bin after refinement when their neighbours sit in it and the bin's markers do not repeat. auto walks down from `--min-contig-size` in bands, each holding as many contigs as are already in play. A contig joins when more of its neighbours sit in the bin than the band's share bar, which starts at one half and rises until fewer than half of the joining marker contigs above it repeat the bin's markers in place, pooled from the cutoff down. The walk stops at the first band where no bar clears, or whose markers could not add a genome. Each bin refuses its joins on its own markers. A band whose marker contigs would join less sequence than `--min-bin-size` is not searched. A length attaches from there whatever the markers say. Shorter contigs take no part |
| `--min-bin-size` | 200000 | Clusters smaller than this are not written as a bin |
| `-t, --threads` | 10 | Threads for rosella and for everything it calls |
| `--no-refine` | off | Skip the splitter, which cuts up the chimeric bins the first clustering leaves behind |
| `--assembly-graph` | none | A GFA whose links join the neighbour graph as extra edges |

Everything else is in `rosella recover --full-help`, grouped by what it touches. The flags under
"Reports" change no bin and only write a table about the run.
