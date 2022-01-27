---
title: Usage
---

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
rosella bin -r scaffolds.fasta --coverage-values coverm.cov -o rosella_bin/ -t 24
```

### Option 2: Use rosella to get coverage values

If you have yet to run `CoverM` then rosella can generate the coverage values for you!
This is especially useful if you have both long and short reads as they can be passed
to rosella in tandem. Don't be afraid to pass multiple samples at once to rosella either,
it can handle it and keep everything in order.

```
rosella bin -r scaffolds.fasta -1 short_s[12345].1.fastq.gz -2 short_s[12345].2.fastq.gz --longreads nanopore.fastq.gz -o rosella_bins/
```

You can keep the BAM files that are created by rosella by including the `--bam-file-cache-directory`
flag. Once the coverage values are calculated, they will be stored in the output direcotry along with
the kmer frequency file.

### Option 3: Use previous results

If you have previously run rosella and something happened like a crash but the coverage values have
already been calculated and stored in your provided output directory then you are safe.
Just run rosella by giving it your output directory and watch as it sorts everything out for you.

```
rosella bin -r scaffolds.fasta -o rosella_bins/ -t 24
```

Rosella checks for previous results by default, but this behaviour can be overridden if you'd like fresh coverage
values by passing the `--force` flag.

## Outputs

The main output for rosella will be a set of MAGs denoted `rosella_bin_X.fna`. How many bins you get depends on your 
samples. Additionally, the kmer frequency table will be present: `rosella_kmer_table.tsv`. And the coverage values if
they were calculated by rosella: `rosella_coverage_values.tsv`. Finally, you'll get a pretty UMAP projection plot coloured
by potential MAG clusters. This plot isn't necessary but it helps you interpret how well rosella partitioned out your contigs.
If you see only a couple of big noisy clusters then maybe something went wrong and you'll want to fiddle with a few of
the UMAP parameters. This is unlikely though, but if you feel like you do need to then please feel free to raise an issue
on this GitHub and I'll answer your question and add my respone to the FAQ to help other users.