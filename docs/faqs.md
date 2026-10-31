# FAQs

**How many samples do I need?**

One works. More is much better. Coverage over a single sample is one number per contig, so two
organisms sitting at the same depth cannot be told apart on it, and composition has to do the
whole job. Every extra sample makes the coverage view harder to fake.

**Do I have to run CoverM first?**

No. Pass reads and rosella maps them itself. strobealign and minimap2 are built in.

**Can I use a depth table from somewhere else?**

If it is in CoverM contig format, yes. That is what `coverm contig -m metabat` writes and what
most pipelines already have lying around. A table holds no window depth though, so pass the BAMs
it was made from when you still have them.

**Why are two runs of the same data identical, and why is there no `--seed` that changes that?**

Because the partition ensemble and the neighbour search take fixed seeds rather than the
command line one. Run variance was measured to sit almost entirely in the partition seed, so it
was fixed rather than left to the user. `--seed`, `--knn-seed` and `--partition-seed` are there
to probe that, not to be tuned.

**What does `--no-refine` turn off?**

The splitter, which cuts up the chimeric bins the first clustering leaves behind. It runs by
default at one round. Rounds two onward measured identical to one.

Leave it on unless you are short of time. The other finishing stages are not affected. Everything named in `--stage-order`, including
the rescue pool and the shed, runs either way.

**Do I need a GPU?**

No. `--gpu` makes training and the big neighbour searches faster, and a run without a GPU does the
same work on the CPU. If `--gpu` finds no GPU, the log says why on a line starting `GPU: none`.

**What are `rosella_bin_unbinned.fna` and `rosella_bin_small_unbinned.fna`?**

Contigs with no bin. The small one holds contigs under `--attach-floor`, which never took part in
binning at all. Between the two of them and the real bins,
every contig in the assembly is written exactly once.

**Can I trust `quality.tsv`?**

For steering a run, yes. The `gtdb_` columns are what the binner decided with, and the `checkm_`
columns, filled with `--checkm`, read like CheckM1 on the same bins. For a paper, score the bins with CheckM2 as well. Markers cannot see contamination
measured in base pairs, so a bin can read clean and still carry foreign sequence.

**Is there a `--markers`, a `--composition-metric`, a `--join-whole`?**

No. Flags that won a comparison were folded into the behaviour and deleted. If you found one in
an old note or an old page, check `rosella recover --full-help` before believing it.

**It says my rosella version is `-dirty`. Does that matter?**

It means the binary was built from a tree with uncommitted changes, so the commit in the version
string does not fully describe it. It runs fine. Do not attribute a result to that commit.
