## UniCoracle: Automated Hierarchical Feature Selection

Open access web server: [micportal.org](https://micportal.org)

Source code: [![DOI](https://zenodo.org/badge/676863744.svg)](https://doi.org/10.5281/zenodo.19050205)

Supplementary Data: [![DOI](https://zenodo.org/badge/1098234753.svg)](https://doi.org/10.5281/zenodo.19050197)


Combines [UniCorP bottom up propagation](https://doi.org/10.1093/ismeco/ycaf174) with [Coracle](https://doi.org/10.1093/bioinformatics/btad749) top down skimming (TDS).
Exploits taxonomic structure of microbiome data.
* UniCorP bottom up approach enriches higher taxonomic levels
* TDS Coracle approach returns to lowest hierarchical level


## Input Format

Requires three input datasets:

### 1. Feature data file (`x`)
A sample-by-feature matrix, where rows represent individual samples and columns represent features (e.g., ASVs or gene IDs). This table should contain raw counts or relative abundances depending on the intended transformation. It must not contain any non-numeric metadata columns.

### 2. Target variable data file (`y`)
A single-column table containing the continuous target variable (e.g., pH, temperature, biomass) for each sample. The index must match the sample IDs in the feature table.

### 3. Hierarchical structure data file (`tax`)
A feature-by-level matrix describing the biological or functional hierarchy of the features. Each row corresponds to a feature (matching the columns in the feature table), and each column represents a hierarchical level, from broad (e.g., Phylum) to fine-grained (e.g., Genus or ASV).
Missing values can be left blank or replaced with the feature ID (recommended) itself. The order of levels should ideally go from highest (left) to lowest (right), although UniCorP can infer and adjust for reversed hierarchies.

### Hierarchy Requirements
* The hierarchy must be a strict tree: each child node should have exactly one parent, and the number of nodes should decrease toward broader levels (e.g., Phylum > Class > Order).
* UniCoracle cannot operate on DAGs (directed acyclic graphs) such as those found in Gene Ontology or KEGG, as these violate the strict parent-child structure.
* Flat hierarchies (where each node has only one child per level) provide no exploitable structure for propagation and are therefore not suitable.

### Format Summary

| File | Orientation | Index | Columns | Notes |
|---|---|---|---|---|
| `x` (Feature) | Samples × Features | Sample IDs | Feature IDs | Numeric only |
| `y` (Target) | Samples × 1 | Sample IDs | Target name | Continuous variable |
| `tax` (Hierarchy)| Features × Levels | Feature IDs | Hierarchical Levels | Strict tree |

### Core Parameters

* `n_features`: Maximum features per hierarchical level - minimum 2 features (to keep compositionality), recommended **100**
* UniCor metric selection: `uc_top_k` or `uc_threshold`, recommended **`uc_top_k` = 100** (like `n_features`)
* Correlation method (`uc_method`): Pearson or **Spearman**
* Transformation (`uc_transformation`): **Relative abundance**, centered log ratio, or raw counts
