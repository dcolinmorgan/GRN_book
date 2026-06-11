# PANDA: Passing Attributes between Networks for Data Assimilation

## Overview

PANDA (Passing Attributes between Networks for Data Assimilation) is a message-passing algorithm for estimating bipartite gene regulatory networks (GRNs). It integrates three data types:

1. **Gene expression** data (samples × genes)
2. **Transcription factor binding motif (TFBM)** prior network (TF → gene edges)
3. **Protein-protein interaction (PPI)** data (TF–TF edges)

PANDA iteratively refines a seed GRN (from motif data) by passing messages between these data layers until convergence, producing edge weights that represent the strength of evidence for each TF–gene regulatory relationship.

## Method

The algorithm works by:
- Initializing a regulatory network from TF binding motif predictions
- Computing a gene co-expression network from expression data
- Iteratively updating edge weights using message passing between the motif network, co-expression network, and PPI network
- Converging on a final network where edge weights reflect integrated evidence across all three data types

The update step uses Tanimoto similarity to measure concordance between data sources and balances "availability" (evidence that a TF regulates a gene) with "responsibility" (evidence from TF protein-complex partners).

## Citation

Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing Messages between Biological Networks to Refine Predicted Interactions. *PLoS ONE*. 2013;8(5):e64832. doi:[10.1371/journal.pone.0064832](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832)

## Software

This chapter uses [netZooPy](https://github.com/netZoo/netZooPy), the Python implementation from the Network Zoo project.

- Documentation: https://netzoopy.readthedocs.io
- PyPI: `pip install netZooPy`
- CLI: `netzoopy panda --e expression.txt --m motif.txt --p ppi.txt --o output.txt`

## Input Data

| File | Format | Description |
|------|--------|-------------|
| Expression | Tab-separated, genes as rows, samples as columns (no header by default) | Gene expression matrix |
| Motif | Tab-separated: `TF  Gene  Weight(0/1)` | Prior regulatory network from motif scans |
| PPI | Tab-separated: `TF1  TF2  Weight` | Protein–protein interactions between TFs |

Toy data for testing is bundled with netZooPy at `tests/ToyData/`.

## Requirements

```
netZooPy>=0.10
pandas
numpy
matplotlib
```

## Links

- [netZoo project](https://netzoo.github.io/)
- [GRAND database](https://grand.networkmedicine.org/) — pre-computed PANDA networks for GTEx tissues
- [Netbooks](https://netbooks.networkmedicine.org/) — Jupyter tutorials from the netZoo team
