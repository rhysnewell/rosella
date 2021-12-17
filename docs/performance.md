---
title: Performance
---

Performance
========

*The following is  breakdown of some preliminary rosella results. This is in now way suggesting you should
use rosella over other tools, instead it is meant to make you consider using rosella in conjunction with other tools.
Each binner has strengths and weaknesses, and when used together via DASTool it is their joint strengths that begin
to shine rather than their weaknesses.*

Testing and benchmarking for rosella is still underway, but the initial results look pretty good.
Comparing performance on the CAMI one low, medium, and high complexity datasets shows that rosella compares as well as,
if not better than other binning tools. In the high complexity dataset this even includes DASTool which is a combination
of the results of all the other binners except for VAMB and rosella.

The results of each benchmark is posted below. Each bin set was filtered down to only include bins that had >=50%
Completeness and <=10% Contamination based on the results of `checkm lineage_wf` v1.1.2.
The bins were than ranked by their completeness and contamination values separately and the results plotted.

### CAMI Low

![CAMI Low](/results/cami_low.png)

Here Rosella manages to outperform every other independent binning algorithm.

### CAMI Medium

![CAMI Med](/results/cami_med.png)

Here rosella generates more bins of higher quality than
the other single binning techniques. DASTool still results in the best results
but the inclusion of rosella in the DASTool algorithm would only be beneficial.

### CAMI High

![CAMI High](/results/cami_high.png)

Here we rosella really take flight, outperforming even DASTool.


### CAMI metagenomes breakdown

| Community | Samples  |   | Inputs |   |
|---|---|---|---|---|
|   | Short read | Long read | Total Gbp | Genomes | Circular elements |
|   CAMI Low	| 1	| 0	| 15 | 40 | 20   |
|   CAMI Med	| 4	| 0	| 40 | 132 | 100 |
|   CAMI High	| 5	| 0	| 75 | 596 | 478 |