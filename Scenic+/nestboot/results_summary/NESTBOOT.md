# NestBoot tests

Tests use all 50 pairs (runs 0–49). Every null `eRegulons.tsv` is from the current eGRN pass; 19 early nulls were rebuilt from their saved `tf_to_gene.tsv` on 22 Oct 2026. The older tables remain beside them as `eRegulons.legacy.tsv` and are not read.

A feature is scored once per run, present or absent. Run i real is paired with run i null. The exact McNemar test asks whether, among runs where the two disagree, the feature is more often present in real. p-values are one-sided. Benjamini–Hochberg FDR is computed separately for TF–gene links and for eRegulon names. Absence from every tested run is not a test.

## TF–gene links

Tested **170,689** links that appear in at least one current real or null run. **35,047** have McNemar q ≤ 0.05 (28,775 at q ≤ 0.01). All 35,047 of the q ≤ 0.05 calls are in the real-greater direction.

Of the q ≤ 0.05 links, **15,413** are also present in at least 25/50 real runs. Those are the links counted per TF in `tf_significant_targets.tsv` (34 TFs). Median real support among these stable discoveries is 40/50; median null support is 0/50.

| TF | significant target genes (q ≤ 0.05 and real support ≥ half) |
|---|---:|
| MEF2C | 2847 |
| ZNF148 | 1749 |
| SOX5 | 1419 |
| BCL11A | 1217 |
| TFDP2 | 966 |
| NFIB | 879 |
| RARB | 868 |
| TCF12 | 798 |
| NFIA | 669 |
| SOX6 | 440 |
| FLI1 | 415 |
| ZNF385D | 411 |
| IKZF1 | 375 |
| ETS2 | 316 |
| NF1 | 297 |

## eRegulons

Tested **891** eRegulon names. **66** have McNemar q ≤ 0.05 (58 at q ≤ 0.01). The unpaired Fisher exact test, same direction, calls **68** at q ≤ 0.05. Agreement on the McNemar discoveries: 65 of 66 also have Fisher q ≤ 0.05.

eRegulons with McNemar q ≤ 0.05, ordered by real inclusion then mean target-gene count:

| eRegulon | real runs | null runs | mean genes real | mean genes null | q McNemar | q Fisher |
|---|---:|---:|---:|---:|---:|---:|
| MEF2C_direct_+/+ | 50 | 8 | 2504.4 | 37.0 | 9.65e-12 | 8.46e-19 |
| MEF2C_direct_+/- | 50 | 5 | 1250.6 | 25.8 | 1.81e-12 | 2.19e-21 |
| BCL11A_direct_+/+ | 50 | 8 | 1221.0 | 48.1 | 9.65e-12 | 8.46e-19 |
| NFIB_direct_+/+ | 50 | 6 | 1026.6 | 22.2 | 2.98e-12 | 1.69e-20 |
| SOX5_direct_+/- | 50 | 5 | 1019.1 | 27.8 | 1.81e-12 | 2.19e-21 |
| NFIA_direct_+/+ | 50 | 9 | 726.3 | 18.1 | 1.76e-11 | 4.82e-18 |
| TCF12_direct_+/+ | 50 | 2 | 667.7 | 18.5 | 3.52e-13 | 1.30e-24 |
| TCF12_direct_+/- | 50 | 1 | 453.7 | 10.0 | 2.64e-13 | 7.51e-26 |
| FLI1_direct_+/+ | 50 | 15 | 413.5 | 36.3 | 9.60e-10 | 6.54e-14 |
| ETV6_direct_+/+ | 50 | 16 | 400.0 | 18.1 | 1.67e-09 | 2.44e-13 |
| IKZF1_direct_+/+ | 50 | 6 | 369.8 | 29.5 | 2.98e-12 | 1.69e-20 |
| SOX5_direct_-/+ | 50 | 0 | 331.4 |  | 2.64e-13 | 4.42e-27 |
| SOX5_direct_+/+ | 50 | 2 | 266.6 | 20.0 | 3.52e-13 | 1.30e-24 |
| RFX4_direct_+/+ | 50 | 5 | 260.2 | 27.2 | 1.81e-12 | 2.19e-21 |
| SOX6_direct_+/+ | 50 | 4 | 236.2 | 20.2 | 1.15e-12 | 2.79e-22 |
| CREB5_direct_+/+ | 50 | 2 | 147.2 | 11.0 | 3.52e-13 | 1.30e-24 |
| RBFOX2_direct_+/+ | 50 | 0 | 84.5 |  | 2.64e-13 | 4.42e-27 |
| ETS2_direct_+/+ | 49 | 12 | 339.1 | 27.9 | 2.03e-09 | 2.32e-14 |
| SOX6_direct_-/+ | 49 | 0 | 226.4 |  | 2.64e-13 | 7.51e-26 |
| RUNX1_direct_+/+ | 49 | 0 | 115.4 |  | 2.64e-13 | 7.51e-26 |
| ABL1_direct_+/+ | 49 | 0 | 13.1 |  | 2.64e-13 | 7.51e-26 |
| ZNF148_direct_+/- | 48 | 2 | 1472.2 | 10.5 | 1.15e-12 | 1.21e-21 |
| SOX4_direct_+/+ | 47 | 13 | 186.9 | 23.8 | 1.67e-09 | 1.44e-11 |
| ZNF385D_direct_+/+ | 46 | 3 | 425.4 | 31.0 | 5.63e-12 | 1.95e-18 |
| NF1_direct_+/- | 44 | 0 | 331.7 |  | 2.98e-12 | 1.69e-20 |
| RFX4_direct_-/+ | 42 | 0 | 45.5 |  | 9.65e-12 | 8.46e-19 |
| NCALD_direct_+/+ | 41 | 0 | 38.4 |  | 1.76e-11 | 4.82e-18 |
| GLI3_direct_+/- | 38 | 1 | 60.0 | 11.0 | 2.70e-10 | 2.32e-14 |
| TFDP2_direct_+/- | 36 | 6 | 1457.3 | 26.5 | 8.35e-07 | 1.34e-08 |
| RORA_direct_+/+ | 36 | 0 | 30.4 |  | 5.19e-10 | 1.76e-14 |
| RARB_direct_+/+ | 35 | 0 | 720.7 |  | 9.60e-10 | 6.54e-14 |
| PHF20_direct_+/+ | 35 | 2 | 167.3 | 16.5 | 2.75e-08 | 2.21e-11 |
| PBX1_direct_+/- | 34 | 0 | 204.4 |  | 1.67e-09 | 2.44e-13 |
| RB1_direct_+/+ | 34 | 0 | 66.0 |  | 1.67e-09 | 2.44e-13 |
| ATF6_direct_+/+ | 31 | 1 | 42.1 | 13.0 | 2.51e-08 | 2.66e-10 |
| RARB_direct_+/- | 30 | 4 | 1281.7 | 54.2 | 9.21e-06 | 4.29e-07 |
| RCOR1_direct_+/+ | 30 | 10 | 283.9 | 42.3 | 2.78e-03 | 7.40e-04 |
| ZNF148_direct_+/+ | 29 | 3 | 1726.7 | 30.0 | 2.41e-06 | 2.20e-07 |
| MEF2C_direct_-/- | 29 | 0 | 112.5 |  | 4.61e-08 | 1.38e-10 |
| ZBTB16_direct_+/+ | 29 | 0 | 94.8 |  | 4.61e-08 | 1.38e-10 |
| TAL1_direct_+/+ | 26 | 7 | 172.7 | 45.3 | 4.27e-03 | 8.34e-04 |
| RFX4_direct_+/- | 26 | 2 | 48.4 | 18.5 | 1.36e-06 | 6.55e-07 |
| FOXN2_direct_+/+ | 24 | 0 | 30.2 |  | 1.36e-06 | 3.57e-08 |
| NKX6-2_direct_+/- | 22 | 5 | 87.7 | 16.6 | 3.88e-03 | 1.84e-03 |
| NKX6-2_direct_+/+ | 22 | 8 | 65.0 | 18.9 | 2.98e-02 | 2.83e-02 |
| CREB5_direct_-/+ | 22 | 0 | 52.7 |  | 5.18e-06 | 2.70e-07 |
| CUX2_direct_+/+ | 21 | 2 | 164.7 | 17.0 | 6.26e-04 | 6.78e-05 |
| RARB_direct_-/+ | 20 | 0 | 20.7 |  | 1.98e-05 | 1.82e-06 |
| EMX2_direct_+/+ | 19 | 5 | 305.9 | 29.2 | 2.98e-02 | 1.36e-02 |
| CUX1_direct_+/+ | 19 | 3 | 100.0 | 10.7 | 3.38e-03 | 1.51e-03 |
| ZNF385D_direct_-/+ | 19 | 0 | 61.9 |  | 3.78e-05 | 4.55e-06 |
| PRRX1_direct_-/+ | 19 | 0 | 44.8 |  | 3.78e-05 | 4.55e-06 |
| SHOX_direct_+/+ | 16 | 1 | 51.9 | 24.0 | 2.40e-03 | 7.00e-04 |
| KLF12_direct_+/+ | 15 | 4 | 1200.8 | 41.0 | 4.99e-02 | 6.04e-02 |
| GLI3_direct_-/+ | 15 | 0 | 25.7 |  | 5.91e-04 | 1.68e-04 |
| MEF2C_direct_-/+ | 14 | 1 | 84.4 | 13.0 | 2.18e-03 | 3.03e-03 |
| GLIS3_direct_-/+ | 14 | 0 | 36.0 |  | 1.13e-03 | 3.94e-04 |
| SOX15_direct_+/+ | 14 | 2 | 33.9 | 20.5 | 2.96e-02 | 1.35e-02 |
| SOX4_direct_-/+ | 13 | 0 | 36.8 |  | 2.18e-03 | 8.55e-04 |
| CUX1_direct_+/- | 12 | 0 | 184.1 |  | 3.88e-03 | 1.84e-03 |
| UBE2K_direct_+/- | 12 | 0 | 35.0 |  | 3.88e-03 | 1.84e-03 |
| ZNF91_direct_+/- | 11 | 0 | 244.0 |  | 7.50e-03 | 4.05e-03 |
| ETV1_direct_-/- | 10 | 0 | 27.8 |  | 1.43e-02 | 8.67e-03 |
| BCL11A_direct_-/+ | 10 | 0 | 17.1 |  | 1.43e-02 | 8.67e-03 |
| EMX2_direct_-/+ | 10 | 0 | 15.5 |  | 1.43e-02 | 8.67e-03 |
| RFX4_direct_-/- | 9 | 0 | 59.6 |  | 2.81e-02 | 1.83e-02 |

## Files

- `results_summary/tf_gene_nestboot.tsv` — every tested TF–gene link, with discordant counts, McNemar p and q.
- `results_summary/eregulon_nestboot.tsv` — every tested eRegulon, McNemar and Fisher.
- `results_summary/tf_significant_targets.tsv` — per-TF count of stable real-enriched target genes.

Inclusion ignores repeated regions for the same TF–gene pair inside one run. Target-gene means for an eRegulon are averaged only over runs in which that eRegulon is present.
