## UniCoracle: Automated Hierarchical Feature Selection

Open access web server: [micportal.org](https://micportal.org)
Source code: [SebastianStaab/UniCoracle](https://github.com/SebastianStaab/UniCoracle)

Combines [UniCorP bottom up propagation](https://doi.org/10.1093/ismeco/ycaf174) with [Coracle](https://doi.org/10.1093/bioinformatics/btad749) top down skimming (TDS).
Exploits taxonomic structure of microbiome data.
* UniCorP bottom up approach enriches higher taxonomic levels
* TDS Coracle approach returns to lowest hierarchical level


## Input Format

Requires three input datasets:

### 1. Feature data file (`x`)
* Sample by feature matrix
* Rows: individual samples
* Columns: features
* Raw counts or relative abundances
* No non numeric metadata columns

### 2. Target variable data file (`y`)
* Continuous variable measurements
* Matches sample set

### 3. Hierarchical structure data file (`tax`)
* Strict tree like feature hierarchy

### Format Summary

| File | Orientation | Index | Columns | Notes |
|---|---|---|---|---|
| `x` (Feature) | Samples × Features | Sample IDs | Feature IDs | Numeric only |
| `y` (Target) | Samples × 1 | Sample IDs | Target name | Continuous variable |
| `tax` (Hierarchy)| Levels × Features | Levels | Feature IDs | Strict tree |

### Core Parameters

* `n_features`: Maximum features per hierarchical level - minimum 2 features (to keep compositionality), recommended **100**
* UniCor metric selection: `uc_top_k` or `uc_threshold`, recommended **`uc_top_k` = 100** (like `n_features`)
* Correlation method (`uc_method`): Pearson or **Spearman**
* Transformation (`uc_transformation`): **Relative abundance**, centered log ratio, or raw counts
