# Gene Regulatory Network Inference: Methods and Protocols

Companion code repository for the forthcoming volume in the Springer [*Methods in Molecular Biology*](https://link.springer.com/series/7651) series.

This repository contains data, notebooks, and reproducible demos for each chapter/method covered in the book. Each folder corresponds to a GRN inference method and includes the materials needed to follow along with the chapter protocols.

## Repository Structure

| Folder | Method | Category |
|--------|--------|----------|
| `Beeline/` | [BEELINE](https://github.com/Murali-group/Beeline) | Benchmarking |
| `BiGSM/` | BiGSM | ML / Deep Learning |
| `Bootstrap/` | NestBoot | Statistical |
| `D-Spin/` | D-Spin | Statistical Mechanics |
| `DeepRIG/` | DeepRIG | ML / Deep Learning |
| `GeneSPIDER/` | GeneSPIDER | Benchmarking / Simulation |
| `GLUE/` | GLUE | Multi-omics Integration |
| `LINGER/` | LINGER | Multi-omics Integration |
| `LiPLike/` | LiPLike | Likelihood-based |
| `merlin/` | MERLIN | Network Inference |
| `RegDiffusion/` | RegDiffusion | Diffusion Models |
| `Scenic+/` | SCENIC+ | Multi-omics (scRNA + scATAC) |
| `scMTNI/` | scMTNI | Multi-omics / Temporal |

## Getting Started

Each method folder contains (or will contain) self-contained demo notebooks and any required input data. To reproduce a chapter's analysis:

1. Clone this repository:
   ```bash
   git clone https://github.com/dcolinmorgan/GRN_book.git
   cd GRN_book
   ```
2. Navigate to the method of interest, e.g.:
   ```bash
   cd Scenic+/region_sets
   jupyter notebook demo.ipynb
   ```
3. Follow the instructions in the chapter-specific `readme.md` or notebook for environment setup and dependencies.

## About the Book

This volume is part of the [*Methods in Molecular Biology*](https://link.springer.com/series/7651) series (Springer), which provides step-by-step protocols for biological research. Our book focuses on computational methods for **Gene Regulatory Network (GRN) inference** from single-cell and bulk transcriptomic data, covering:

- Perturbation-based approaches
- Machine learning and deep learning methods
- Statistical and likelihood-based inference
- Multi-omics integration (scRNA-seq + scATAC-seq)
- Benchmarking frameworks

## Contributing

Each chapter author maintains their own method folder. If you find issues reproducing a demo, please open an issue referencing the specific method/chapter.

## License

Please refer to individual chapter folders for licensing of code and data. The book content is published by Springer Nature.
