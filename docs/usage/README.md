Getting started
========

Every flag is documented in the binary itself. `rosella recover -h` lists the ones in
everyday use and `rosella recover --help` lists all of them, the same for `rosella refine`
and `rosella score`.

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
Optionally, you can also pass a bin quality table with `--bin-quality` to limit refinement to only MAGs contaminated above `--split-contamination`. A CheckM1, CheckM2 or AMBER table is read as-is.

```bash
rosella refine -r scaffolds.fasta -d metabat_bins/ -x fna -C coverm.cov -o refined_bins/ -t 24
```
OR
```bash
rosella refine -r scaffolds.fasta -d metabat_bins/ -x fna -1 short_s[12345].1.fastq.gz -2 short_s[12345].2.fastq.gz --longreads nanopore.fastq.gz -o refined_bins/ -t 24
```

Every bin lands in `--output-directory`, named `rosella_<--bin-tag>_N.fna` and numbered from
zero. Bins the refiner did not change are written out alongside the ones it split, under new
numbers, so the output names do not carry over from the input. The directory has to be empty
of `.fna` files: a second run into a used one is refused rather than appended to.

## Outputs

The main output for rosella will be a set of MAGs denoted `rosella_bin_X.fna`. How many bins you get depends on your 
samples. Alongside them are the per-stage timings as `timings.tsv` and a marker-based quality table as
`quality.tsv`. And the coverage values if they were calculated by rosella: `coverage.tsv`.
`--write-kmer-table` also keeps the composition table, named for the k it counted, as in
`kmer_frequencies.k4.tsv`. It runs to hundreds of megabytes on a large assembly, so it is off by
default and a later run over the same assembly reuses whatever it finds.
If you see only a couple of big clusters then something went wrong. `--partition` and
`--n-neighbours` are the two knobs that change the shape of the answer most; raise an issue
on this GitHub and I'll answer and add the response to the FAQ.
