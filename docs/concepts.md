# How it works

Rosella bins by cutting a graph, not by clustering points. The sections below follow the order
the stages run in.

## Features

Two views of every contig. Composition is its tetranucleotide frequencies. Depth is its mean
coverage in each sample. When rosella reads BAMs or reads, depth also holds every 100 bp window
along the contig.

Depth is the stronger of the two and it gets weaker the fewer samples you have. Over a single
sample it is one number per contig, so two organisms at the same depth are indistinguishable on
it and composition has to carry the separation alone. That is the single largest thing you
control from outside: more samples is more signal.

The single copy markers are annotated alongside, from the genes rosella calls on every contig
the partition takes.

## The learned space

A network maps each contig's composition and depth to a point. The neighbour graph, the
partition, the rescue pool and attach read contigs through those points.

With window depth the network learns without labels, while the markers are still being
annotated. Each contig has eight fixed views. One is the whole contig. Seven are random stretches
of it, each at least 1 kb long. Each training step shows the whole contig and `--stretches` of
its stretches. A stretch reads its depth from the 100 bp windows it covers. Half the views also
borrow a random share of the depth of a contig with close composition, as reads from a sister
strain would. The space still has to place such a view with its own contig, so it learns not to
part contigs on exact depth alone. Training pulls the views of one contig together and pushes
them away from every other contig's views. Depth has a tower of its own whose output joins the
composition, so a contig's depth pattern weighs the same whatever its abundance. Training runs
for up to `--epochs` passes over the contigs. It stops early once `--patience` passes bring no
real fall in loss.

A coverage table holds no windows, so each contig has only one view. The run first bins on the
raw k-mers and depth. It runs one pass at a starting weight, fits the coverage weight from that
pass's near complete bins, and runs a second pass at that weight. The near complete bins of the
second pass become labels for a smaller network. The network's output then stands in for the
k-mers, and depth joins it at a fitted weight. This recovers fewer genomes than window depth
does.

`--gpu` trains the network and embeds the contigs on the GPU. `--space-cache` keeps a space
trained on window depth for a later run over the same inputs.

A few late checks still read the raw tables. These are the splitter's trigger, the replicon and
passenger tests, and the fold of marker-free bins.

## The neighbour graph

Each contig is joined to its 100 nearest neighbours in the learned space. The neighbourhoods are
folded into a fuzzy simplicial set, which turns a hard neighbour list into weighted edges.
`--n-neighbours` sets the width. On the CPU the big searches are approximate and seeded to a
fixed value, so they are reproducible without being exhaustive. Under `--gpu` they are exact.

With `--assembly-graph` the links from a GFA join that graph as extra edges, weighted by
`--assembly-graph-weight`.

## The partition, which is an ensemble

The graph is cut several ways rather than once. Leiden runs across a ladder of resolutions,
label propagation runs alongside it, and both run again at each of `--partition-seeds` seeds.
No one of those cuts is the answer.

## Arbitration

Each Leiden ladder hands over the two rungs with the most bins over the bar. Neighbouring rungs
often tie within a few bins, and each keeps whole some genomes the other cuts, so one pick would
be a coin flip. Label propagation makes one cut per seed and hands that over. Every community in
those cuts becomes a candidate, and the bins are chosen one at a time. Single copy markers score
each candidate, the best ones claim their contigs, and a candidate whose contigs are already
taken is offered what is left. When a larger candidate holds at least half of the winner and
scores level with it inside one standard error of the markers, the larger one is taken. A
genome's contigs without a marker look the same to the markers as a gap.

This is why the markers matter so much. They are the only thing in the pipeline that knows what
a genome is, as against what a dense part of a graph is.

`rosella ensemble` adds each other binner's bins to these candidates as one more rung, so a bin
from any of them can win its genome. Each winner then goes to the binners' vote. A contig joins it
when at least half of the binners whose good bins overlap it hold that contig, and leaves when
most binners keep it in a good bin elsewhere. The markers cannot see a genome's contigs that carry
none, and the other binners can. The vote runs once more on the finished bins before attach.

## Settling the coverage weight

Depth joins the learned space at a weight, and the partition runs more than once. This holds
with window depth too, where the space already carries depth through its own tower. Each pass
refits the weight from the near complete bins of the pass before. A new pass is kept only when
its marker score beats the kept pass by more than the spread between passes, the one standard
error rule (Breiman et al. 1984). The walk stops at the first pass that fails, at a
weight it has already tried, or after eight passes. The log names these passes `pass 1` to
`pass 8`.

The fit compares each contig's two halves. Each half is embedded in the learned space at its own
depth, drawn at half the contig's span.

## Markers that run in families

The bins of the kept pass that are at least half complete show which markers this community
carries in more than one copy. A second genome doubles every marker it carries at once. A marker
doubled far more often than the rest is a family with copies of its own, so from then on it
counts once.

## Genome scale

The near complete bins of the kept pass measure the smallest genome in the run: a bin's bases
over its completeness. With no near complete bin the run takes 400 kb. Every size the cycle cares
about follows from that one number. A bin claiming a share of a genome must hold at least that
share of the smallest genome's bases, so a dense cluster of marker genes cannot pass for a
genome. A cut may not leave a piece smaller than the smallest genome at the lowest completeness a
bin is reported at, and a contig on its own is written as a bin from that length. None of this is
the output floor. `--min-bin-size` only stops clusters under 10 kb being written, so a small
genome or phage the partition holds together is kept.

## The outlier pool

Long contigs that no bin claimed get one partition of their own. Every community it finds becomes
a bin, for the stages below to judge.

## The splitter

The splitter measures each bin's spread in depth and in the learned space. A bin whose spread
sits above most bins' is clustered again. So is a bin over `--max-bin-size`. A cut stands when
the pieces are tighter than the whole and two of them are large enough to be part of a genome.
The splitter also peels off a long contig that sits far from the rest of its bin, unless the
markers say the contig belongs. It runs one round, because rounds two onward measured identical
to one. `--no-refine` skips it.

## The refine cycle

The cycle runs its stages in the order `--stage-order` names, by default:

- **dissolve** is the rescue pool. Bins that miss the completeness and contamination bars go
  back into a pot with the unbinned contigs. The pot's neighbour graph is cut from the run's
  graph and partitioned again. Each of the `--dissolve-rounds` rounds halves the neighbour width.
  A genome buried inside a dense graph then gets a chance to form its own community. Candidates that clear the strict
  bar are adopted first. Looser bars then take what is left, down to half complete. What a claim
  leaves of a bin is kept as a bin, and a bin the pool broke and then made nothing better of is
  put back the way it was.
- **join** offers the scorer whole pairs of bins to fuse. It is off by default and `--join`
  turns it on. The markers cannot tell two halves of one genome from two genomes, so about half
  of what it fuses on real data is two genomes.
- **shed** unbins a contig whose every marker the bin already holds. It only thins bins that
  miss the quality bars, since thinning a bin that is already good only costs it sequence. It is
  skipped when at least half of the bins handed to the pool already clear the bars.

## Attach

Contigs from `--attach-floor` up to `--min-contig-size` took no part in the partition. They now
join the bin their neighbours sit in. Each bin weighs the marker contigs about to join it. If
their repeated markers outweigh the markers they fill, the bin refuses every join.

## Publishing

Before the bins are written:

- Replicons and passengers leave their bins. [Output](output.md) says what each is.
- Bins under `--min-bin-size` are dropped, and their contigs go to the leftovers.
- Bins that code less densely than any prokaryote in the run are set apart as eukaryotes.
- A bin under 5 per cent complete folds into a good bin whose depth across the samples it
  matches. A good bin here is at least 50 per cent complete and at most 10 per cent contaminated.

Every contig is then written exactly once, to a bin or to one of the leftover files.

## What it cannot do

Three things are measured and worth knowing before you go looking for a knob.

A long contig placed in the wrong community by the partition is very hard to recover
afterwards. No stage that moves contigs between finished bins gets much of it back, because
neither composition nor single sample coverage separates a passenger from the genome it is
riding inside.

Single copy markers reach only a small share of contigs. A contig carrying no marker cannot be
judged by them at all, whatever the stage.

The marker table cannot see contamination measured in base pairs. A bin can read 99 per cent
complete and 2 per cent contaminated and still be carrying a meaningful fraction of foreign
sequence, because that sequence carried no duplicated marker.
