# How it works

Rosella bins by cutting a graph, not by clustering points. The stages below run in this order.

## Features

Two views of every contig. Composition is the tetranucleotide frequency table, and coverage is
whatever the depth table holds, one column per sample.

Coverage is the stronger of the two and it gets weaker the fewer samples you have. Over a single
sample it is one number per contig, so two organisms at the same depth are indistinguishable on
it and composition has to carry the separation alone. That is the single largest thing you
control from outside: more samples is more signal.

## The neighbour graph

Each contig is joined to its 100 nearest neighbours under a distance that combines the two
views, and the neighbourhoods are folded into a fuzzy simplicial set, which is what turns a hard
neighbour list into weighted edges. `--n-neighbours` sets the width. The neighbour search is
approximate and seeded to a fixed value, so it is reproducible without being exhaustive.

With `--assembly-graph` the links from a GFA join that graph as extra edges, weighted by
`--assembly-graph-weight`.

## The partition, which is an ensemble

The graph is cut several ways rather than once. Leiden runs across a whole ladder of
resolutions, label propagation runs alongside it, and the ladder is built again at each of
`--partition-seeds` seeds. No one of those cuts is the answer.

## Arbitration

Every community any rung proposed becomes a candidate, and the bins are chosen one at a time
rather than by picking a winning rung. Single copy markers score each candidate, the best ones
claim their contigs, and a candidate whose contigs are already taken is offered what is left.

This is why the markers matter so much. They are the only thing in the pipeline that knows what
a genome is, as against what a dense part of a graph is.

## The rescue pool

Bins that miss the completeness and contamination bars are not kept as they are. They go back
into a pot with the unbinned contigs, that pot is embedded and partitioned again, and candidates
that clear the bars are adopted. The search narrows over `--dissolve-rounds` neighbour widths, so
a genome buried inside a dense graph gets a chance to form its own community.

A bin the pool broke and then made nothing better of is put back the way it was.

## Finishing

The splitter runs first, cutting up the bins the markers call chimeric. One round, because
rounds two onward measured identical to one. `--no-refine` skips it.

The refine cycle then runs its stages in the order `--stage-order` names, by default:

- **dissolve**, the rescue pool above
- **join**, offering the scorer whole pairs of bins to fuse, including an incomplete partner for
  a bin the markers already call whole
- **recruit**, off by default, letting a bin short of the bars take single contigs off a neighbour
- **audit**, reading each short contig against its own neighbourhood
- **shed**, unbinning a contig whose every marker the bin already holds, and only from bins that
  miss the quality bars, since thinning a bin that is already good only costs it sequence

Then the bins are written.

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
