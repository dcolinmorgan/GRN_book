# Seed consensus on the SCENIC+ runs

GRNBoost2, the TF-to-gene step of SCENIC+, returns the same network for a fixed expression matrix and a fixed seed, and a different network when the seed changes. These 50 measured runs use one multiome and seeds `42 + run`. Nothing else changes: cells are not resampled, and region-to-gene links, cistromes, and the ranking database are fixed. The 50 shuffled runs permute each gene’s expression across cells before GRNBoost2, with seeds `42 + run + 10000`. They are the control.

A link’s **seed support** is the fraction of the 50 runs that contain it. A link that the data supports should come back under most seeds. A link that only reflects one random draw of the solver should not.

## Self-links and sparsity

Self-links are removed first. GRNBoost2 ranks a TF predicting itself first in every run, measured and shuffled, because the two expression columns are the same vector, and the gene shuffle does not break that edge. Those 700 edges are not regulatory links. They were the spike at support 1.0 in earlier weight-ranked plots.

Each run is then cut to a few links per node by GRNBoost2 importance. Two directed readings are reported: the 3 strongest regulators of each gene (87,495 edges per run) and the 3 strongest targets of each TF (2,100 edges per run). Degrees 4 and 5 and the uncut table are in the last section.

## The measured networks recur; the shuffled ones do not

At 3 regulators per gene, two measured runs share a Jaccard overlap of 0.30 before any cutoff. Two shuffled runs share 0.003. 62% of measured links appear in at least 2 seeds, against 13% of shuffled links. Restricting to links that recur raises the measured overlap: 0.58 for links in at least half the seeds, 0.80 for links in at least 40 of 50. No shuffled link appears in 8 or more seeds. 924 measured links are in all 50.

At 3 targets per TF the measured runs agree less (overlap 0.13, against 0.0001 for the shuffle), and the same pattern holds: half of the measured links appear in at least 2 seeds, against under 1% of shuffled links. Overlap is 0.50 at support 0.5 and 0.79 at support 0.8. Six links are in all 50 seeds. The top targets of a TF move a lot between seeds, but far less than chance.

## Cutoff from the shuffled control

The cutoff is set by the shuffled runs, not by a chosen percentage of seeds. For a support value B,

FDR(B) = (shuffled links with support > B) / (measured links with support > B)

and the cutoff is the smallest B at which this is at most 0.05. This is Equation 1 of Morgan et al. (2019), used here as an empirical false-discovery rate against the shuffled control.

### 3 regulators per gene

The ratio first falls to 0.05 at **B = 0.06**: the link is present in **at least 4 of 50** seeds.

| Support greater than | Measured links | Shuffled links | FDR |
|---:|---:|---:|---:|
| 0 (seen at least once) | 551,565 | 3,828,642 | 6.94 |
| 2/50 | 270,227 | 48,828 | 0.181 |
| **3/50 (the cutoff)** | **230,264** | **4,069** | **0.018** |
| 4/50 | 203,680 | 288 | 0.0014 |
| 6/50 | 168,861 | 1 | 6×10⁻⁶ |
| 7/50 | 156,133 | 0 | 0 |
| 40/50 | 19,494 | 0 | 0 |
| 50/50 | 924 | 0 | 0 |

### 3 targets per TF

The ratio first falls to 0.05 at **B = 0.02**: the link is present in **at least 2 of 50** seeds.

| Support greater than | Measured links | Shuffled links | FDR |
|---:|---:|---:|---:|
| 0 (seen at least once) | 27,570 | 104,309 | 3.78 |
| **1/50 (the cutoff)** | **13,775** | **683** | **0.050** |
| 2/50 | 9,525 | 8 | 0.0008 |
| 3/50 | 7,319 | 0 | 0 |
| 40/50 | 80 | 0 | 0 |
| 50/50 | 6 | 0 | 0 |

The cutoff is low because a shuffled matrix almost never produces the same edge twice. The set it keeps is a union across seeds, larger than any one run: 230,264 links against 87,495 per run at 3 regulators per gene, and 13,775 against 2,100 at 3 targets per TF. Report seed support with every link. Then a reader can separate the links that only pass the control from the core that recurs in most seeds.

## 3–5 links per node, and every edge

![3 and 5 links per node](plots/fig2_degree_3_and_5.png)

![All TF–gene edges](plots/fig2_all_edges.png)

The 2×2 is degrees 3 and 5. The top row is regulators per gene, the bottom row targets per TF. The second figure is every TF–gene edge GRNBoost2 wrote, with no degree cut. The x-axis is seed support. Each point is the raw fraction of links at that support (1/50, 2/50, …). The curves are not smoothed; the line only joins those points. The y-axis runs from 0 to 0.2, and the one-seed peak sits above that frame. Blue is measured, red is shuffled, and the black line is that panel’s cutoff. `fig2_grid.py` writes the values to `fig2_degree_and_full.tsv`.

| Graph | Edges per run | Cutoff B | Measured links above it | Shuffled links above it | FDR | Shuffled frequency is 0 from |
|---|---:|---:|---:|---:|---:|---:|
| 3 TFs per gene | 87,495 | 0.06 | 230,264 | 4,069 | 0.018 | 0.16 |
| 4 TFs per gene | 116,660 | 0.06 | 296,611 | 11,917 | 0.040 | 0.16 |
| 5 TFs per gene | 145,825 | 0.08 | 320,337 | 3,391 | 0.011 | 0.18 |
| 3 targets per TF | 2,100 | 0.02 | 13,775 | 683 | 0.050 | 0.08 |
| 4 targets per TF | 2,800 | 0.04 | 12,326 | 12 | 0.001 | 0.08 |
| 5 targets per TF | 3,500 | 0.04 | 15,111 | 19 | 0.001 | 0.08 |
| All edges | ~2.4 million | 0.42 | 1,566,575 | 55,314 | 0.035 | 0.78 |

On the degree-cut panels the shuffled frequency drops to zero by support 0.08–0.18, while the measured curve continues to 1.0. The uncut table is the one where shuffled links still recur: shuffled frequency stays positive until support 0.78, and the cutoff moves out to 22 of 50 seeds. Keeping every edge is what lets shuffled links recur. Cutting to a few links per node removes that tail.

## Why this is not NestBoot

NestBoot (Morgan et al., *Bioinformatics* 35:1026–1032, 2019) estimates support from bootstrap resamples of independent experiments. It repeats that inner loop in an outer loop with a fresh experiment-wise shuffle, and keeps a link only when it clears the Equation 1 cutoff in every outer run. That last rule is why the FDR in the paper’s Figure 6 is lower than the rate Equation 1 asks for. The method assumes the inference is fixed once the input is fixed, so all variation comes from the data.

These runs do not meet those assumptions. There is one multiome, so there are no experiments to resample. The measured matrix is the same in all 50 runs, and only the GRNBoost2 seed changes, so the support measures the solver’s variation, not the sample’s. There is one set of runs and no outer loop. Earlier tables split the 50 seeds into 5 groups of 10 or 10 groups of 5 (`inner_bootstrap_fdr.tsv`). That only repartitions the same seeds. The shuffle permutes each gene within the one dataset rather than shuffling experiments.

A NestBoot on this dataset would draw cells with replacement, keep the GRNBoost2 seed fixed, and repeat Q inner and R outer times on measured and shuffled data. One GRNBoost2 fit took about 1.5 hours here, so even R = 5 and Q = 20 means 200 fits. That was not run.

Other tables on disk answer different questions and are not part of this result: a sweep of GRNBoost2 rank cutoffs (`sparsity_fdr.tsv`), the eRegulon-table support (`scenicplus_support_fdr.txt`), and a per-link paired test (`tf_gene_nestboot.tsv`).
