# Usage

## Let rosella map the reads

Pass reads instead and rosella maps short reads with strobealign 0.19 and long reads with
`minimap2 -x map-ont`, both built in. `--mapper minimap2-sr` maps short reads with minimap2
instead. `--longread-mapper` picks another long read preset: `minimap2-pb`, `minimap2-hifi`,
`minimap2-lr-hq` or `minimap2-no-preset`. Long and short
reads can go in together, and so can as many samples as you have. More samples is better: coverage over one sample is a
single number per contig, and a lot of what the binner can separate comes from how that number
moves across samples.

```bash
rosella recover -r assembly.fasta \
  -1 sample_{1,2,3}.1.fq.gz -2 sample_{1,2,3}.2.fq.gz \
  --longreads nanopore.fq.gz \
  -o rosella_bins/ -t 24
```

Each sample's BAM is kept in `bams/` under the output directory, so a later run can take them
back with `-b`. It stays in the mapper's order, with every hit that scores as well as the best.
A kept BAM older than the assembly or its reads is mapped again, and so is one whose `@PG` line
names another mapper or other `--minimap2-parameters`. Read files must be plain or
gzipped. Each gzipped read file is read to its end alongside the mapping. A damaged one stops the run
and keeps no BAM, since the mapper would quietly stop at the damage and drop every read after it.
An index larger than a quarter of the free memory is built in parts that each fit, sized from
a slice of the assembly, and minimap2 merges the hits across parts. A later sample reuses any
index already split at least as finely. Other `--minimap2-parameters` get an index of their own,
since minimap2 keeps the k-mer settings an index was built with. Mapping then runs at full speed instead of swapping.
strobealign builds its index for each sample in seconds and keeps none. An index larger than
two thirds of the free memory is split into parts of about equal bases. Every part maps every
read, and each mate keeps its best hit from any part as its primary, beside every hit level with it.

## BAMs you already have

BAMs skip the mapping, short and long reads alike. They can be coordinate sorted, such as the
ones `coverm make` or aviary keeps, or in the mapper's order as rosella keeps them:

```bash
rosella recover -r assembly.fasta -b binning_bams/*.bam -o rosella_bins/ -t 24
```

Rosella reads each BAM once. Every counted alignment covers its reference span, with no
identity filter and no end trimming. A coordinate sorted BAM counts each primary alignment
whole. A BAM in the mapper's order keeps a read's equally good hits together, so a read whose
best score is shared by several contigs is split between them by their depth, refitted until
the split settles: a strain at 10x takes ten times the shared reads of one at 1x. Both mates of a
pair are scored together. Each mate keeps only its hits on the contigs where the two mates' best
scores add up highest. A pair whose mates share no contig keeps its hits as they are. That one
pass gives the coverage table and the depth of every 100 bp window. The learned space reads each
view of a contig from those windows. A coverage table alone carries no windows, so the space
learns instead from the near complete bins of a first pass on the raw k-mers and depth. That
recovers fewer genomes. Give BAMs or reads when you have them.

The depth read from each BAM is kept in `depth/` under the output directory, or wherever
`--depth-directory` points:

| File | What it holds |
|---|---|
| `contigs.tsv` | Each contig and its length, in the BAMs' header order |
| `<sample>.depth.zst` | One zstd stream, little-endian, in `contigs.tsv` order: every contig's mean depth as f64, then every variance as f64, then every primary alignment count as u32, then the mean depth of every 100 bp window as half floats, contig by contig. A contig's last window covers what is left of it |
| `coverage.tsv.gz` | Every sample's mean and variance in CoverM's contig table layout, written by `rosella make` |
| `<sample>.summary.tsv` | Reads in the BAM, how many mapped, how many tied across contigs, mean read length and mapped bases. Unmapped reads count only when the BAM keeps them, as rosella's own BAMs do |
| `<sample>.ties.zst` | Contig pairs that share tied reads, as little-endian u32 triples: first contig, second contig, reads, in `contigs.tsv` order |
| `<sample>.links.zst` | Contig pairs joined by one read's supplementary alignment, in the same layout. Long reads that span two contigs land here |
| `<sample>.bam.key` | A key from the length, head and tail of the BAM the depth was read from |

The read count says how sure a mean is. Neighbouring bases share reads, so the variance across a
contig cannot say that, and a long read contig at low depth can rest on a few dozen reads.

A later run given the same BAMs and the same `--depth-directory` reads the kept depth and skips
the BAM, for every sample whose BAM still matches its key. A BAM copied in with `rsync -a` keeps
its old mtime, so the key reads its bytes instead. Other tools can read the depth too,
so binners run beside rosella need not walk the BAMs again.

## Depth without keeping the BAMs

`rosella make` writes the depth directory and nothing else, so the BAMs can be deleted as soon as
it finishes. It can run once per sample, each run adding to the directory and rewriting
`coverage.tsv.gz` from every sample in it. After runs in parallel, one more with no reads or BAMs
rebuilds that table alone. `recover` reads every sample in the directory in name order:

```bash
rosella make -b binning_bams/*.bam -o depth/ -t 24
rm binning_bams/*.bam
rosella recover -r assembly.fasta --depth-directory depth/ -o rosella_bins/ -t 24
```

Reads work too, with `-r` for the assembly to map them to. A minimap2 index of the assembly is
kept beside the BAMs, so a run per sample indexes it once. strobealign indexes again each run.

The composition table is only kept with `--write-kmer-table`, since it runs to hundreds of
megabytes on a large assembly and costs seconds to rebuild.

## A coverage table

A CoverM contig table, made with `coverm contig -m metabat`, works when the BAMs are gone:

```bash
rosella recover -r assembly.fasta -C coverage.tsv -o rosella_bins/ -t 24
```

A table gives each contig one depth per sample and no windows. The learned space then trains on
the near complete bins of a first pass. That recovers fewer genomes than BAMs or reads do.

## The GPU

`--gpu` trains the learned space, embeds the contigs and runs the big neighbour searches on the
GPU. Under it those searches are exact rather than approximate, so a run reads a few bins
differently from the same run on the CPU. A run that finds no GPU logs why and carries on on the
CPU. [Installation](installation.md#the-gpu) says what each platform needs.

## Refining bins that already exist

`rosella refine` takes bins from any binner and runs them through recover's refine cycle. The
bins stand in for recover's first partition. Everything after it runs as recover runs it: the
markers, the learned space, the splitter, the rescue pool, the shed, attach and publishing. It
takes every flag recover takes except `--partition-seeds` and the probes of the first partition.
It also takes the bins:

```bash
rosella refine -r assembly.fasta -d metabat_bins/ -x fa -C coverage.tsv -o refined_bins/ -t 24
```

The run holds to the contigs in the bins. A contig in no bin is never read or written, so refining a
few bins leaves the rest of the assembly alone. A bin the refine cycle leaves whole is written
whatever its size and is never folded into another, since the input binner already chose to keep it.
`--min-bin-size` only floors the pieces refine cuts. Short contigs in the bins wait for attach, as
in recover. A contig named in two bins stays in the first. Without `--marker-cache`, genes are
called on the bins' contigs alone. With it, the annotation of the whole assembly is read, or written
once for recover to share. With `--space-cache`, the space trained over the whole assembly is read
in place of one trained on the bins' contigs, and trained and kept first if no run has kept it.

Pass a bin quality table with `--bin-quality` to aim the splitter: a CheckM1, CheckM2 or AMBER
table is read as it is, and bins over `--split-contamination`, 15 per cent by default, are always
split candidates.

`--output-directory` has to hold no `.fna` files already. A second run into a used directory is
refused rather than appended to. Bins the refiner did not change are written out beside the ones
it split, renumbered, so the output names do not carry over from the input.

## Binning with several binners' bins

`rosella ensemble` runs recover with the bins of other binners offered beside its own partition.
Each entry is one binner's output, gzipped or not. A directory holds one `.fa`, `.fna` or
`.fasta` file per bin. A file is a table of contig and bin, tab separated, with lines starting `@`
or `#` skipped, so a CAMI binning or a contig-to-bin table reads without writing the bins out:

```bash
rosella ensemble -r assembly.fasta --depth-directory depth/ -o ensemble_bins/ \
    -d metabat_bins/ semibin_bins/ quickbin.tsv -t 24
```

Each binner's bins are one more arm of the per-bin combination that picks recover's first
partition. A bin is taken where its markers score higher than every bin it shares contigs with, so
a binner can win one genome and lose the next. At an equal score, the bin more binners offer whole
goes first. When a larger bin holds at least half of the winner and reads level with it inside one standard
error of the markers, and at least as many binners offer it whole, the larger bin is taken instead,
since a genome's own contigs without a marker look the same to the markers as a gap. The bins on
offer differ by whole genomes cut into short contigs, so here a gene cut by a contig end counts as
a second copy wherever it covers the same part of its marker as another copy. The two halves of
one split gene still count once. The winner's contigs
are then put to a vote among the binners whose good bins share most of it. A contig joins when at
least half of them hold it, and leaves when most binners keep it in a good bin elsewhere.
Rosella's own rungs vote as one binner. Everything after that runs as recover runs it, over the
whole assembly, except that the refine cycle's bins are voted on once more against every
binner's bins before attach. The output is recover's. A
contig named in two bins of one binner stays in the first. Across binners it is a contest. A bin
with no markers cannot win a contest, so a virus another binner kept alone is left to the replicon
and circle paths. Contigs under `--min-contig-size` wait for attach like any other short contig.
There each binner that bins one votes, beside its neighbours, for the bin its own bin maps to. A
binner that leaves a short contig out has no say on it, since binners bin few short contigs.

`ensemble_info.tsv` says, per binner, how many bins and bp it offered and how many written bins
trace back to it.
`bin_sources.tsv` traces each written bin. A cached learned space
and marker annotation serve `ensemble` as they serve `recover`, so a run after `recover` with the
same `--space-cache` and `--marker-cache` trains nothing and calls no genes.

## Scoring bins

`rosella score` calls genes on every contig in the bins it is given, searches them for the single
copy markers and writes a completeness and contamination table. No assembly, no gold standard
and no reference database:

```bash
rosella score -d rosella_bins/ -x fna -o quality.tsv -t 24
```

`--proteins <file>` also writes every called gene's protein, named `<contig>_<gene number>` as
prodigal names them. Scoring the whole assembly this way hands a bin refiner the same gene calls
rosella used:

```bash
rosella score -f assembly.fasta.gz -o quality.tsv --proteins proteins.faa -t 24
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

## Reusing the learned space

Training the learned space is most of what is left once the markers are cached. `--space-cache
<dir>` keeps the trained space, keyed on the assembly's contigs, each sample's depth file, the
training flags and a digest of the code that trains it. A later recover or refine with the same
key reads it and trains nothing, so a rosella that changed only how bins are drawn reuses it:

```bash
rosella recover -r assembly.fasta --depth-directory depth/ -o rosella_bins/ \
    --marker-cache cache/ --space-cache cache/
rosella refine -r assembly.fasta --depth-directory depth/ -d rosella_bins/ -o refined/ \
    --marker-cache cache/ --space-cache cache/
```

When no entry matches, refine trains over every contig recover would read and keeps that space,
so each later refine in a loop reads it. Only a space trained on window depth is kept. A run
from a coverage table trains its space every time.

## Reproducibility

Two runs of one assembly with the same inputs on the same machine give the same bins. The
partition ensemble and the neighbour search take fixed seeds rather than `--seed`, so there is
nothing to choose and no reason to run three times. A space trained on a GPU can, rarely, drift
by a few bins between runs. Pass `--space-cache` when two runs must match exactly. `--seed`, `--knn-seed` and `--partition-seed` exist to probe that, not
to be tuned.

## Flags worth knowing

| Flag | Default | What it does |
|---|---|---|
| `--min-contig-size` | auto | Contigs at or above this are partitioned. auto reads 1500. A shorter cutoff takes nothing attach does not, and costs wall |
| `--attach-floor` | auto | Contigs from this length up to `--min-contig-size` join a bin after refinement when their neighbours sit in it and the bin's markers do not repeat. auto walks down from `--min-contig-size` in bands, each holding as many contigs as are already in play. A contig joins when more of its neighbours sit in the bin than the band's share bar, which starts at one half and rises until fewer than half of the joining marker contigs above it repeat the bin's markers in place, pooled from the cutoff down. The walk stops at the first band where no bar clears, or whose markers could not add a genome. Each bin refuses its joins on its own markers. A band whose marker contigs would join less sequence than the smallest piece of a genome is not searched. A length attaches from there whatever the markers say. Shorter contigs take no part |
| `--min-bin-size` | auto | Clusters smaller than this are not written as a bin. auto is 10 kb, the least a replicon needs to be written alone, so a small genome the partition holds together is written. Cutting and pooling are held to the genome scale instead |
| `-t, --threads` | 10 | Threads for rosella and for everything it calls |
| `--gpu` | off | Train the learned space and run the big neighbour searches on the GPU when there is one |
| `--fast-train` | off | Train the learned space mostly on contigs it has not learned yet. About a quarter of the training time, for a few bins fewer |
| `--no-refine` | off | Skip the splitter that cuts and peels chimeric bins |
| `--assembly-graph` | none | A GFA whose links join the neighbour graph as extra edges |

Everything else is in `rosella recover --full-help`, grouped by what it touches. The flags under
"Reports" change no bin and only write a table about the run.
