# Reproducibility, and what a single SCENIC+ run lets through

Fifty paired bootstraps. Each real run is GRNBoost2 `TF_to_gene` on the observed expression (seed `42 + run`) followed by eGRN with the fixed `region_to_gene` table from `real_0`. The paired null shuffles each gene’s expression independently (seed `42 + run + 10000`) and repeats that same eGRN. A TF–gene link is counted once per run, however many regions repeat it. All 50 null `eRegulons.tsv` files are from the current eGRN pass (scenicplus 1.0a2, git `840dab85de3044846234c157e327b6ea0290abfb`).

## A single run still emits a large null network

The eGRN thresholds are not a null filter. Every real run beats its paired null (50/50; Wilcoxon signed-rank p = 8.9×10⁻¹⁶), but the null is not empty.

| | real | null |
|---|---:|---:|
| triplets, median (min–max) | 47,111 (37,551–59,705) | 3,995 (685–8,384) |
| TF–gene links, median | 19,165 | 2,022 |
| eRegulon names, mean | 53.7 | 73.7 |

A shuffled run still returns about **2,000 TF–gene links** and **more eRegulon names** than the real run. Those null eRegulons are small (fewer genes and regions per name). Median TF–gene rho is 0.288 in real runs and 0.076 in null runs, so the links that survive a shuffle are weaker, but they do pass the same pipeline. Median fold in triplet count is 12.3×. Paired real-vs-null Jaccard is 0.001 (mean overlap **16** links, range 0–42). One null replicate therefore cannot be subtracted from one real network link by link: the passing sets barely intersect, and both are large.

## Bootstrap recurrence removes the links a single run would keep

Real TF–gene sets are only moderately stable (mean pairwise Jaccard **0.424**). Null sets are not (Jaccard **0.002**). Of the links present in a typical real run (~19,244):

| recurrence | links in the union | share of one real run | also McNemar q ≤ 0.05 |
|---|---:|---:|---:|
| ≥10/50 | 28,924 | 84% | all of them |
| ≥25/50 | 15,413 | 63% | all of them |
| ≥40/50 | 7,989 | 39% | all of them |
| 50/50 | 2,077 | 11% | all of them |

So a single published eGRN keeps a large fringe that this bootstrap does not reproduce. About **37%** of the TF–gene links in one real run are present in fewer than half of the 50 seeds. About **61%** miss a 40/50 bar. Only **11%** are in every real run. Every link that clears 10/50 real runs is also enriched versus the null at q ≤ 0.05 (35,047 links at q ≤ 0.05 overall; median null support among the ≥25/50 set is 0/50).

The one-seed fringe is the part the null also produces. **19,100** real links appear in exactly one run. None of them reach q ≤ 0.05. On the null side, recurrence collapses: 95,360 distinct null links appear at least once, **142** reach 5/50, **60** reach 10/50, **5** reach 25/50, and **none** reach 40/50. A link that the shuffle can regenerate does not come back in the next shuffle. Requiring the same link across seeds is what removes it. A single real-versus-null comparison of network size would not, because both sides pass the eGRN cutoffs.

The same split shows up at eRegulon names. 891 names appear at least once. **17** are in every real run and **27** are in at least 40/50; all 27 of those have q ≤ 0.05. **No** null eRegulon reaches 40/50, even though a null run carries more names (mean 74 vs 54). 86 null names reach 10/50, and only 5 of those are still called real-enriched. Names that merely clear one eGRN run, including most null names, do not survive the bootstrap.

## What to treat as reproduced

For this benchmark, a TF–gene link is reproduced when it is present in at least 25 of 50 real runs. That set has 15,413 links, all with q ≤ 0.05, median real support 40/50 and median null support 0/50. The TFs that contribute the most of these targets are MEF2C (2,847), ZNF148 (1,749), SOX5 (1,419), BCL11A (1,217), and TFDP2 (966). An eRegulon name is reproduced on the same rule when its real inclusion is high and its null inclusion is not; the 66 names with McNemar q ≤ 0.05 are listed in `NESTBOOT.md`. Fisher exact tests agree on 65 of those 66.

These filters do not add an expression or motif threshold beyond SCENIC+. They count how often the same stochastic GRNBoost2 step rediscovers the link when region-to-gene is held fixed, and how often the gene-wise shuffle does too.
