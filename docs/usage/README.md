Getting started
========

## Inputs

Rosella can be run multiple different ways in order to make using it as easy as possible.

### Option 1: Precomputed coverage values

The most likely situation is that you are running rosella in conjunction with other
binning algorithms, like metabat2. If so, then you likely already have coverage values precomputed
using `coverm contig` using `-m metabat`.

To perform mag recovery:
```
rosella recover -r scaffolds.fasta --coverage-file coverm.cov -o rosella_bin/ -t 24
```

### Option 2: Use rosella to get coverage values

If you have yet to run `CoverM` then rosella can generate the coverage values for you!
This is especially useful if you have both long and short reads as they can be passed
to rosella in tandem. Don't be afraid to pass multiple samples at once to rosella either,
it can handle it and keep everything in order.

```
rosella recover -r scaffolds.fasta -1 short_s[12345].1.fastq.gz -2 short_s[12345].2.fastq.gz --longreads nanopore.fastq.gz -o rosella_bins/
```

Rosella does not keep the BAM files it maps into. Once the coverage values are calculated they
are stored in the output directory as `coverage.tsv`, alongside the kmer frequency file, and
either can be passed back in with `--coverage-file` or `--kmer-frequency-file` to skip that
work on a later run.

### Option 3: Refining the results of other binners

Rosella can also be used to refine the results of other binning algorithms. The required input for this process is:
    - The original assembly FASTA file
    - MAGs from Rosella or another binning algorithm
    - Coverage values for the original assembly OR a set of reads to calculate them with
Optionally, you can also provide the results of CheckMv1 or CheckMv2 (or AMBER) to limit refinement to only MAGs that are contaminated above the max contmiantion threshold.

```bash
rosella refine -r scaffolds.fasta -d metabat_bins/ -x fna -C coverm.cov -o refined_bins/ -t 24
```
OR
```bash
rosella refine -r scaffolds.fasta -d metabat_bins/ -x fna -1 short_s[12345].1.fastq.gz -2 short_s[12345].2.fastq.gz --longreads nanopore.fastq.gz -o refined_bins/ -t 24
```

The output of this process will be two folders, `refined_bins/` and `unchanged_bins/`. The former will contain the refined MAGs and the latter will contain the MAGs that were not refined (either because they were below the contamination threshold or no change occurred after refinement).

## Outputs

The main output for rosella will be a set of MAGs denoted `rosella_bin_X.fna`. How many bins you get depends on your 
samples. Additionally, the kmer frequency table will be present: `kmer_frequencies.tsv`. And the coverage values if
they were calculated by rosella: `coverage.tsv`.
If you see only a couple of big noisy clusters then maybe something went wrong and you'll want to fiddle with a few of
the UMAP parameters. This is unlikely though, but if you feel like you do need to then please feel free to raise an issue
on this GitHub and I'll answer your question and add my respone to the FAQ to help other users.
