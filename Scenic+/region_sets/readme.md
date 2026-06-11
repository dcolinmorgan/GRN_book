# DeepWiki Q&A with Code Context for Repository: aertslab/pySCENIC
## Q1
how can we use this to create a GRN using function calls not docker or snakemake
### Answer
The user wants to know how to create a Gene Regulatory Network (GRN) using pySCENIC through Python function calls (not Docker or Snakemake). Let me explore the repository to understand the API and workflow.

Let me get more details about the key functions - particularly `modules_from_adjacencies`, `prune2df`, `df2regulons`, and `aucell` - to give a complete picture of the Python API workflow.

I have enough information to provide a comprehensive answer about how to use pySCENIC to create a GRN using Python function calls. Let me compile the complete workflow.


Here is the complete Python API workflow to build a GRN with pySCENIC — no Docker or Snakemake needed.

---

## Required inputs

1. Expression matrix (cells x genes) — CSV, loom, or h5ad
2. TF list — a `.txt` file with one TF per line
3. cisTarget ranking databases — `.feather` files (download from [cistargetDBs](https://resources.aertslab.org/cistarget/))
4. Motif annotation table — a `.tbl` file

---

## Step 0: Imports and setup

```python
import os, glob, pickle
import pandas as pd
from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2          # or genie3

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
``` 
pySCENIC:46-63 

---

## Phase I: GRN inference (co-expression modules)

```python
# Load expression matrix: rows=cells, columns=genes
ex_matrix = pd.read_csv("expression.tsv", sep='\t', header=0, index_col=0).T

# Load TF names
tf_names = load_tf_names("tf_list.txt")

# Run GRNBoost2 to get TF-target adjacencies
adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

# Convert adjacencies to co-expression modules (candidate regulons)
modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
```

`grnboost2` returns a DataFrame with columns `TF`, `target`, `importance`. pySCENIC:131-162 

`modules_from_adjacencies` applies several thresholding strategies (percentile cutoffs, top-N targets, top-N regulators) and uses Pearson correlation to split activating vs. repressing modules. pySCENIC:267-299 

---

## Phase II: Prune with cisTarget (motif enrichment → regulons)

```python
# Load ranking databases
db_fnames = glob.glob("path/to/databases/*.genes_vs_motifs.rankings.feather")
dbs = [RankingDatabase(fname=f, name=os.path.splitext(os.path.basename(f))[0])
       for f in db_fnames]

# Prune modules: keep only targets with cis-regulatory motif support
with ProgressBar():
    df = prune2df(dbs, modules, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")

# Convert enriched motif table to Regulon objects
regulons = df2regulons(df)

# Optionally save
df.to_csv("motifs.csv")
with open("regulons.pkl", "wb") as f:
    pickle.dump(regulons, f)
``` 
pySCENIC:168-180 

---

## Phase III: AUCell — score regulon activity per cell

```python
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
# auc_mtx: DataFrame of shape (n_cells x n_regulons)
auc_mtx.to_csv("auc_matrix.csv")
``` 
pySCENIC:206-209 

---

## Important notes

- **If running as a script** (not a Jupyter notebook), wrap your code in `if __name__ == '__main__':` to avoid multiprocessing issues. pySCENIC:39-41 

- **Parallelism**: `grnboost2` uses Dask internally. For local multiprocessing without Dask, you can use the `arboreto_with_multiprocessing` script directly via Python's `multiprocessing.Pool`. pySCENIC:127-188 

- **Sparse matrices**: Pass a sparse CSC matrix to `grnboost2` along with `gene_names` separately if memory is a concern. pySCENIC:96-104 

- The `prune2df` step is the most resource-intensive. You can pass `client_or_address` to use a Dask distributed cluster if needed. pySCENIC:184-188

## pySCENIC Python API Workflow for GRN Creation
This codemap traces the complete pySCENIC workflow for creating Gene Regulatory Networks using Python function calls, covering the four main phases: GRN inference with GRNBoost2/GENIE3 [1d], module creation from adjacencies [2a], cisTarget-based pruning with motif enrichment [3b], and AUCell regulon activity scoring [4a].
### 1. GRN Inference Phase
Traces the initial network inference using arboreto's GRNBoost2 or GENIE3 algorithms to generate TF-target adjacencies from expression data
### 1a. Import arboreto functions (`tutorial.rst:56`)
Import load_tf_names and grnboost2 from the arboreto package for GRN inference
```text
from arboreto.utils import load_tf_names
    from arboreto.algo import grnboost2
```
### 1b. Load expression matrix (`tutorial.rst:85`)
Load single-cell expression matrix with cells as rows and genes as columns
```text
ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T
```
### 1c. Load transcription factors (`tutorial.rst:97`)
Load the list of known transcription factors from a text file
```text
tf_names = load_tf_names(MM_TFS_FNAME)
```
### 1d. Run GRN inference (`tutorial.rst:133`)
Execute GRNBoost2 algorithm to infer TF-target weighted adjacencies from expression data
```text
adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)
```
### 1e. CLI GRN inference call (`pyscenic.py:106`)
CLI implementation showing the actual method call with optional Dask client for parallelization
```text
network = method(
                expression_data=ex_mtx,
                tf_names=tf_names,
                verbose=True,
                client_or_address=client,
                seed=args.seed,
            )
```
### 2. Module Creation Phase
Traces conversion of TF-target adjacencies into co-expression modules using thresholding strategies and correlation-based activation/repression classification
### 2a. Create modules from adjacencies (`tutorial.rst:162`)
Convert adjacencies DataFrame into a list of Regulon objects representing co-expression modules
```text
modules = list(modules_from_adjacencies(adjacencies, ex_matrix))
```
### 2b. Module creation function signature (`utils.py:267`)
Function definition showing parameters for thresholding strategies and minimum gene requirements
```text
def modules_from_adjacencies(
    adjacencies: pd.DataFrame,
    ex_mtx: pd.DataFrame,
    thresholds=(0.75, 0.90),
    top_n_targets=(50,),
    top_n_regulators=(5, 10, 50),
```
### 2c. Calculate TF-target correlations (`utils.py:365`)
Add Pearson correlation column to classify targets as activating or repressing based on expression patterns
```text
adjacencies = add_correlation(
                adjacencies,
                ex_mtx,
                rho_threshold=rho_threshold,
                mask_dropouts=rho_mask_dropouts,
            )
```
### 2d. Filter by activation type (`utils.py:376`)
Separate activating from repressing modules and optionally keep only activating ones
```text
activating_modules = adjacencies[adjacencies[COLUMN_NAME_REGULATION] > 0.0]
        if keep_only_activating:
            modules_iter = iter_modules(
                activating_modules, frozenset([ACTIVATING_MODULE])
            )
```
### 2e. Filter and add TF to modules (`utils.py:399`)
Filter modules by minimum gene count, add the TF itself to each module, and return final list
```text
return list(filter(lambda m: len(m) >= min_genes, map(add_tf, modules_iter)))
```
### 3. cisTarget Pruning Phase
Traces motif enrichment analysis using cisTarget ranking databases to prune targets and create final regulons with cis-regulatory evidence
### 3a. Load cisTarget databases (`tutorial.rst:104`)
Load feather-format ranking databases containing genome-wide motif rankings
```text
db_fnames = glob.glob(DATABASES_GLOB)
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
```
### 3b. Run motif enrichment pruning (`tutorial.rst:171`)
Execute distributed motif enrichment analysis to identify modules with cis-regulatory footprint support
```text
with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)
```
### 3c. Set up transformation function (`prune.py:413`)
Configure the function that converts modules to enriched motif dataframes using AUC-first implementation
```text
transformation_func = partial(
        modules2df,
        module2features_func=module2features_func,
        weighted_recovery=weighted_recovery,
    )
```
### 3d. Execute distributed calculation (`prune.py:424`)
Run the pruning calculation using either custom multiprocessing or Dask framework
```text
return _distributed_calc(
        rnkdbs,
        modules,
        motif_annotations_fname,
        transformation_func,
        aggregation_func,
```
### 3e. Convert to regulons (`tutorial.rst:175`)
Transform the enriched motifs dataframe into final Regulon objects with scored target genes
```text
regulons = df2regulons(df)
```
### 4. AUCell Scoring Phase
Traces calculation of regulon activity scores per cell using Area Under the recovery Curve (AUC) analysis
### 4a. Calculate AUCell scores (`tutorial.rst:208`)
Compute regulon enrichment AUC values for each cell using multiprocessing
```text
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)
```
### 4b. AUCell function signature (`aucell.py:185`)
Main function for calculating regulon activity with configurable AUC threshold and worker count
```text
def aucell(
    exp_mtx: pd.DataFrame,
    signatures: Sequence[Type[GeneSignature]],
    auc_threshold: float = 0.05,
```
### 4c. Create gene rankings (`aucell.py:46`)
Convert expression matrix to genome-wide rankings per cell by shuffling and ranking genes by expression
```text
return (
        ex_mtx.sample(frac=1.0, replace=False, axis=1, random_state=seed)
        .rank(axis=1, ascending=False, method="first", na_option="bottom")
        .astype(DTYPE)
        - 1
    )
```
### 4d. Execute AUC calculation (`aucell.py:206`)
Create rankings and pass to aucell4r for the actual enrichment computation with optional multiprocessing
```text
return aucell4r(
        create_rankings(exp_mtx, seed),
        signatures,
        auc_threshold,
        noweights,
        normalize,
        num_workers,
    )
```
