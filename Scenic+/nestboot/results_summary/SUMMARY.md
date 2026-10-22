# NestBoot SCENIC+ summary

Fifty paired bootstraps. Each pair re-runs stochastic GRNBoost2 `TF_to_gene` (seed `42 + run` for real, `42 + run + 10000` for null) and rebuilds eGRNs with the fixed `region_to_gene` table from `real_0`. The null shuffles expression independently within each gene before GRNBoost2. An eRegulon triplet is one row of `eRegulons.tsv` (region, gene, TF).

## Runs

All **50** pairs are on the current eGRN pass. Nineteen early null tables were rebuilt from saved `tf_to_gene.tsv` files on 22 Oct 2026; those previous outputs are `eRegulons.legacy.tsv` and are not used here.

## Size: real versus null

On the 50 current pairs, real triplets exceed the paired null in **50/50** runs (Wilcoxon signed-rank, greater, statistic 1275, p = 8.882e-16).

| | real triplets | null triplets | fold (real/null) | real TF–gene links | null TF–gene links |
|---|---:|---:|---:|---:|---:|
| mean | 47,889.0 | 3,883.4 | 18.5 | 19,243.9 | 1,992.0 |
| sd | 4,917.4 | 1,947.7 | 16.0 | 1,938.0 | 942.3 |
| median | 47,111.0 | 3,994.5 | 12.3 | 19,165.0 | 2,022.0 |
| min | 37,551.0 | 685.0 | 4.9 | 14,805.0 | 365.0 |
| max | 59,705.0 | 8,384.0 | 67.0 | 23,952.0 | 3,953.0 |

| | real mean (min–max) | null mean (min–max) |
|---|---:|---:|
| TFs | 40.0 (32–49) | 64.6 (12–130) |
| genes | 9054.9 (8117–9952) | 1835.3 (360–3488) |
| regions | 30778.1 (25968–36289) | 3633.1 (676–7575) |
| eRegulons | 53.7 (41–67) | 73.7 (12–155) |

No extended triplets appear in the current runs (mean extended fraction 0.000 real, 0.000 null). Positive `regulation` fraction averages 0.942 real and 0.983 null. Median TF-to-gene importance is 1.167 real vs 1.806 null; median rho is 0.288 real vs 0.076 null. Null networks are not just smaller copies of the real ones: they contain more TFs and more eRegulon names, but each covers far fewer genes and regions. Real TF–gene rho is higher; null importance among the triplets that remain is not lower.

## Stability of TF–gene links

Mean pairwise Jaccard of TF–gene sets across the 50 real runs is **0.424** (median 0.423). Across the matched null runs it is **0.002** (median 0.002).

Paired real-vs-null Jaccard on the current runs averages **0.001** (median 0.001), with a mean overlap of 16 links.

Unique TF–gene links seen in at least k of the 50 current runs:

| k | real | null |
|---:|---:|---:|
| 1 | 77,872 | 95,360 |
| 5 | 41,461 | 142 |
| 10 | 28,924 | 60 |
| 25 | 15,413 | 5 |
| 40 | 7,989 | 0 |
| 50 | 2,077 | 0 |

Links in ≥25/50 real runs and ≤2 null runs: **15,396**. Links in ≥40/50 real runs and in no null run: **7,702**.

Union of TF–gene links in the current runs is 170,689 (77,872 real, 95,360 null). `tf_gene_recurrence.tsv` keeps links with real support ≥5 or any null support (135,453 rows). Support columns count only the 50 current runs.

## eRegulon recurrence

891 eRegulon names appear at least once in the current runs. 17 are in every current real run; 0 are in every current null run. 27 real eRegulons reach ≥40/50, against 0 null eRegulons.

Fifteen eRegulons with the highest real support (ties broken by lower null support):

| eRegulon | real runs | null runs |
|---|---:|---:|
| RBFOX2_direct_+/+ | 50 | 0 |
| SOX5_direct_-/+ | 50 | 0 |
| TCF12_direct_+/- | 50 | 1 |
| CREB5_direct_+/+ | 50 | 2 |
| SOX5_direct_+/+ | 50 | 2 |
| TCF12_direct_+/+ | 50 | 2 |
| SOX6_direct_+/+ | 50 | 4 |
| MEF2C_direct_+/- | 50 | 5 |
| RFX4_direct_+/+ | 50 | 5 |
| SOX5_direct_+/- | 50 | 5 |
| IKZF1_direct_+/+ | 50 | 6 |
| NFIB_direct_+/+ | 50 | 6 |
| BCL11A_direct_+/+ | 50 | 8 |
| MEF2C_direct_+/+ | 50 | 8 |
| NFIA_direct_+/+ | 50 | 9 |

## Files

- `results_summary/per_run_counts.tsv` — one row per bootstrap, with `cohort` (`current` or `legacy_null`).
- `results_summary/tf_gene_recurrence.tsv` — TF–gene support in the current runs.
- `results_summary/eregulon_recurrence.tsv` — eRegulon-name support in the current runs.

Counts use presence in a run, not how many regions repeat the same TF–gene pair inside that run.
