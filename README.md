# pySCE

pySCE: Single Cell Entropy Scoring in Python

## Overview

pySCE computes entropy-based metrics for single-cell transcriptomic data. It provides two complementary tools:

1. **Expression Entropy** (`pysce.score`) — Scores each cell's transcriptional entropy using a protein-protein interaction (PPI) network and Markov chain entropy, based on the SCENT framework. GPU-accelerated via PyTorch.

2. **Angular Velocity Entropy** (`pysce.score_angular_velocity_entropy`) — *Experimental.* Scores how coherent or disordered the RNA velocity field is in each cell's local neighborhood. The hypothesis is that this metric may capture dynamic differences that expression entropy alone cannot — for example, tumor-initiating cells (TICs) and normal stem cells can have similar expression entropy but potentially distinct velocity dynamics. This tool is under active development and validation; see [Status](#status) below.

## Installation

```bash
pip install -e .
```

### Requirements

- `scanpy >= 1.8.2`
- `torch >= 2.0.0`
- `scvelo >= 0.2.4`
- `scikit-learn >= 1.0`

## Usage

### Expression Entropy

```python
import pysce

# Score expression entropy (stores in adata.obs['entropy'])
pysce.score(adata)
```

### Angular Velocity Entropy

Requires RNA velocity to be precomputed (e.g., via scVelo) or spliced/unspliced layers for automatic estimation.

```python
import pysce

# Option A: velocity already computed (adata has velocity_umap in .obsm)
pysce.score_angular_velocity_entropy(adata, basis='umap')

# Option B: let pySCE run scVelo for you (needs spliced/unspliced layers)
pysce.ensure_velocity(adata, basis='umap')
pysce.score_angular_velocity_entropy(adata, basis='umap')
```

**Output:**
- `adata.obs['angular_velocity_entropy']` — per-cell score in [0, 1] (normalized), or NaN for cells with zero velocity
- `adata.uns['angular_velocity_entropy_params']` — run metadata including parameters, cell counts, and the normalization denominator

**Interpretation:**
- **Low entropy** (~0): neighbors' velocity vectors are coherent, pointing in similar directions (committed trajectory)
- **High entropy** (~1): neighbors' velocity vectors are scattered (disordered/multipotent neighborhood)

### Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `basis` | `'umap'` | Embedding for neighbor lookup (`'umap'` for 2D, `'pca'` for higher-dim) |
| `n_neighbors` | `30` | Number of nearest neighbors per cell |
| `n_bins` | `8` | Angular histogram bins |
| `normalize` | `True` | Normalize to [0, 1] using dimension-aware max entropy |

### Dimension-Aware Normalization

When `normalize=True`, raw entropy is divided by the theoretical maximum for that dimensionality. For 2D this is `log2(n_bins)` (uniform on circle). For higher dimensions, it accounts for the concentration of measure — random unit vectors in high-dim space have pairwise angles concentrated near pi/2, so the raw max entropy is lower. This makes scores directly comparable across 2D and high-dim bases.

Reference values (n_bins=8):

| Dimensions | Max entropy (bits) | Context |
|---|---|---|
| 2 (UMAP) | 3.00 | Uniform on circle |
| 10 | 1.89 | Moderate concentration |
| 30 (PCA) | 1.21 | Strong concentration near pi/2 |
| 50 | 1.05 | Very tight concentration |

## Status

**Expression Entropy** is based on the established SCENT framework and is stable.

**Angular Velocity Entropy** is speculative and under active validation. The core question is whether velocity-based entropy provides biological insight that expression entropy and existing tools (e.g., scVelo's velocity confidence) do not. Early results on synthetic data confirm the math is sound, and testing on the scVelo pancreas endocrinogenesis dataset shows that the metric captures meaningful structure — actively differentiating populations (Ngn3 high EP, Pre-endocrine) show low angular entropy (coherent flow), while quiescent or terminal populations show higher entropy (disordered neighborhoods). This is consistent with known biology but the interpretation differs from the initial hypothesis: low entropy marks active transitions, not necessarily "committed" fates.

Open questions:
- Does angular velocity entropy distinguish TICs from normal stem cells in practice?
- Does it provide information beyond what scVelo confidence scores already capture? (Preliminary: the two are correlated at r ~ -0.48 but not redundant.)
- How sensitive is the metric to RNA velocity estimation noise, neighbor count, and bin count?

This tool is provided as-is for exploration and hypothesis generation. Interpret results in the context of your specific dataset and biological question.

## License

GNU GPLv3
