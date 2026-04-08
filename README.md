# pySCEv

A fork of [pySCE](https://github.com/dchary/pysce) (Single Cell Entropy Scoring in Python) that adds an experimental RNA velocity entropy module and tests whether it provides useful biological signal — particularly in cases where expression entropy alone cannot distinguish between cell states (e.g., tumor-initiating cells vs. normal stem cells).

The original `pysce.score` expression entropy tool is unchanged. Everything under `_velocity.py` is new.

## What this fork adds

**Angular Velocity Entropy** (`pysce.score_angular_velocity_entropy`) — *Experimental.* For each cell, measures how coherent or disordered the RNA velocity vectors are in its local neighborhood. The hypothesis is that this may capture dynamic differences that expression entropy cannot: two cell populations with identical expression entropy could have very different velocity coherence if one is actively differentiating (coherent flow) while the other is dynamically disorganized.

This is speculative — see [Status](#status) and `examples/` for validation results so far.

## Installation

```bash
pip install -e .
```

### Requirements

- `scanpy >= 1.8.2`
- `torch >= 2.0.0`
- `scvelo >= 0.2.4` (for velocity estimation; optional if you precompute velocity)
- `scikit-learn >= 1.0`

## Usage

### Original Expression Entropy (from pySCE)

```python
import pysce

# Score expression entropy (stores in adata.obs['entropy'])
pysce.score(adata)
```

See the [original pySCE repo](https://github.com/dchary/pysce) for details.

### Angular Velocity Entropy (new in this fork)

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
- **Low entropy** (~0): neighbors' velocity vectors are coherent, pointing in similar directions (active transition / committed trajectory)
- **High entropy** (~1): neighbors' velocity vectors are scattered (disordered / quiescent / multipotent neighborhood)

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

## Examples

See `examples/` for validation figures:

- **`examples/tier1_synthetic/`** — Synthetic data where the ground truth is known. Confirms the math works: aligned velocities score ~0, bifurcation ~0.33, random ~1.0.
- **`examples/tier2_pancreas/`** — scVelo pancreas endocrinogenesis dataset (Bastidas-Ponce et al. 2019). Shows the metric captures meaningful biological structure on real data.

## Status

**Angular Velocity Entropy** is speculative and under active validation. The core question is whether velocity-based entropy provides biological insight that expression entropy and existing tools (e.g., scVelo's velocity confidence) do not.

**What we know so far:**
- Tier 1 (synthetic): Math is sound. Clean separation between aligned (0), bifurcation (0.33), and random (1.0) velocity fields.
- Tier 2 (pancreas): The metric captures meaningful structure — actively differentiating populations (Ngn3 high EP, Pre-endocrine) show low angular entropy (coherent flow), while quiescent or terminal populations show higher entropy (disordered neighborhoods). This is consistent with known biology but the interpretation differs from the initial hypothesis: low entropy marks active transitions, not necessarily "committed" fates. Anti-correlates with scVelo confidence (r ~ -0.48) but is not redundant.

**Open questions:**
- Does angular velocity entropy distinguish TICs from normal stem cells in practice?
- Does it provide information beyond what scVelo confidence scores already capture?
- How sensitive is the metric to RNA velocity estimation noise, neighbor count, and bin count?

This tool is provided as-is for exploration and hypothesis generation. Interpret results in the context of your specific dataset and biological question.

## License

GNU GPLv3
