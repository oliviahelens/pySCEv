# pySCEv

A fork of [pySCE](https://github.com/dchary/pysce) (Single Cell Entropy Scoring in Python) that adds an experimental RNA velocity entropy module and tests whether it provides useful biological signal — particularly in cases where expression entropy alone cannot distinguish between cell states (e.g., tumor-initiating cells vs. normal stem cells).

The original `pysce.score` expression entropy tool is unchanged. Everything under `_velocity.py` is new.

## What this fork adds

**Angular Velocity Entropy** (`pysce.score_angular_velocity_entropy`) — *Experimental.* For each cell, measures how coherent or disordered the RNA velocity vectors are in its local neighborhood. The hypothesis is that this may capture dynamic differences that expression entropy cannot: two cell populations with identical expression entropy could have very different velocity coherence if one is actively differentiating (coherent flow) while the other is dynamically disorganized.

This is speculative — see [Status](#status) and `validation/` for validation results so far.

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

## How It Works

### Expression Entropy (pySCE) — Markov Chain Entropy on a PPI Network

The original pySCE method scores each cell's "potency" by measuring the entropy of a random walk on a protein-protein interaction (PPI) network, weighted by that cell's gene expression.

**Step 1: Build a cell-specific transition matrix.**
Given a PPI adjacency matrix `A` (genes x genes, binary edges) and a cell's expression vector `x` (one value per gene), construct a weighted transition matrix:

```
W_ij = A_ij * x_i * x_j
P_ij = W_ij / sum_j(W_ij)
```

Each entry `P_ij` is the probability of transitioning from gene `i` to gene `j` in one step of a random walk. Genes with higher expression in the cell get more weight, so the random walk preferentially visits highly expressed genes and their PPI neighbors.

**Step 2: Compute the stationary distribution.**
The stationary distribution `pi` (the long-run fraction of time spent at each gene) is computed as:

```
pi_i = x_i * (A * x)_i / sum_j(x_j * (A * x)_j)
```

This weights each gene by its expression times the total expression of its PPI neighbors — genes that are both highly expressed and well-connected in the active network get more weight.

**Step 3: Compute Markov chain entropy.**
The entropy rate of the Markov chain is:

```
S = -sum_i pi_i * sum_j P_ij * log(P_ij)
```

This is the expected Shannon entropy of the next-step distribution, averaged over the stationary distribution. When expression is spread across many PPI-connected genes, the random walk has many possible paths and entropy is high (stem-like / multipotent). When expression is concentrated on a few genes, the walk is constrained and entropy is low (differentiated).

**Step 4: Normalize.**
The score is divided by `log(lambda_1)` where `lambda_1` is the largest eigenvalue of the PPI adjacency matrix. This normalizes the score to [0, 1] relative to the theoretical maximum entropy of the network.

---

### Angular Velocity Entropy (pySCEv) — Shannon Entropy of Local Velocity Directions

The new velocity entropy metric scores each cell by measuring how disordered the RNA velocity directions are in its local neighborhood.

**Step 1: Get velocity vectors and find neighbors.**
For each cell `i`, take its RNA velocity projected into embedding space (e.g., UMAP) as a 2D vector `v_i`. Find the `k` nearest neighbors in embedding space (default `k=30`).

**Step 2: Compute angles.**

*2D embeddings (UMAP):* Convert each neighbor's velocity to an absolute angle:

```
theta_j = atan2(v_j_y, v_j_x),   theta in [0, 2pi)
```

*Higher-dimensional embeddings (PCA):* Compute all pairwise angles among the `k` neighbor velocity vectors using cosine similarity:

```
theta_jl = arccos( (v_j . v_l) / (|v_j| * |v_l|) ),   theta in [0, pi]
```

This gives `k*(k-1)/2` pairwise angles. The pairwise approach is necessary because in high dimensions there is no natural "absolute angle" reference.

**Step 3: Bin and compute Shannon entropy.**
Histogram the angles into `B` equal-width bins (default `B=8`). Compute the Shannon entropy of the resulting distribution:

```
p_b = count_b / sum(counts)
H = -sum_b p_b * log2(p_b)
```

If all neighbors point the same way, one bin gets all the counts and `H = 0`. If directions are uniformly spread, all bins are equal and `H = log2(B)`.

**Step 4: Dimension-aware normalization.**
In 2D, the maximum entropy is simply `log2(B)` (uniform on a circle). In higher dimensions, random unit vectors don't produce uniform pairwise angles — due to *concentration of measure*, angles cluster near `pi/2` as dimensionality grows. The theoretical pairwise angle distribution is:

```
f(theta) proportional to sin^(d-2)(theta),   theta in [0, pi]
```

We integrate this density over each bin to get the expected bin probabilities for random vectors, then compute the Shannon entropy of that distribution. This is the maximum entropy for that dimensionality. Dividing by it normalizes scores to [0, 1] regardless of whether you use UMAP (2D) or PCA (30D).

**Cells with zero velocity** (magnitude < 1e-10) get NaN — entropy of direction is undefined when there is no direction.

## Examples

See `validation/` for figures:

- **`validation/synthetic_validation/`** — Synthetic data where the ground truth is known. Confirms the math works: aligned velocities score ~0, bifurcation ~0.33, random ~1.0.
- **`validation/pancreas_endocrinogenesis/`** — scVelo pancreas endocrinogenesis dataset (Bastidas-Ponce et al. 2019). Shows the metric captures meaningful biological structure on real data.
- **`validation/spatial_chicken_heart/`** — SIRV-imputed Visium (Mantri et al. 2021, day 14). Tests whether physically adjacent neighbors give different answers than UMAP neighbors. They do: r = 0.39 between the two.
- **`validation/parameterization/`** — k neighbors × bin count sweep on pancreas. Median Spearman r = 0.86 across settings; cell-type ordering stable.

## Status

**Angular Velocity Entropy** is speculative and under active validation. The core question is whether velocity-based entropy provides biological insight that expression entropy and existing tools (e.g., scVelo's velocity confidence) do not.

**Notes so far:**
- **Synthetic:** Clean separation between aligned (0), bifurcation (0.33), and random (1.0) velocity fields.
- **Pancreas:** Low entropy marks actively differentiating populations (Ngn3 high EP, Pre-endocrine) — coherent flow, not "committed" fates as the initial hypothesis guessed. Quiescent/terminal populations score higher.
- **Not redundant with existing metrics:** anti-correlates with scVelo confidence (r ~ -0.48, ~77% unique variance), nearly uncorrelated with expression entropy (r ~ 0.02).
- **Correlated with simpler mean angular deviation, r ~ 0.80 on pancreas.** Entropy captures distribution shape; mean angle only captures the first moment. The ~20% unique variance comes from multimodal neighborhoods (genuine bifurcations).
- **Spatial chicken heart (day 14 Visium, 1967 spots):** spatial neighbors and UMAP neighbors give meaningfully different answers (r = 0.39), so physical adjacency captures structure transcriptional adjacency doesn't. In spatial context, entropy and mean deviation diverge more (r = 0.52 vs 0.80) — more cells fall in the multimodal regime where entropy's distribution-shape sensitivity matters. Overlaying Mantri's anatomical regions: lowest entropy in the **left ventricular wall** (trabecular + endocardium, compact myocardium + septum) — the coordinated cardiomyocyte differentiation zones. Highest entropy in **valves and atria** — known developmental mixing zones where multiple divergent trajectories share physical space. Consistent with the pancreas result: low entropy = coordinated differentiation wave, high entropy = multiple divergent trajectories in proximity.
- **Permutation control (N=100) on the chicken heart:** shuffling velocity vectors across spots and recomputing confirms the region ranking is real (Kruskal-Wallis H 329 vs null 66 ± 46, +5.7σ) but reframes the "big smooth tissue patches" claim — about two-thirds of the raw spatial autocorrelation (Moran's I null ≈ 0.47 out of observed 0.74) comes from the k-NN scoring kernel's own smoothing, regardless of biology. The biological signal is the anatomical alignment, not the patchiness itself. See `validation/spatial_chicken_heart/` for details.
- **Parameterization is stable:** median Spearman r = 0.86 across a 5×5 sweep of k and bin count on pancreas; cell-type ordering preserved.

**Open questions:**
- Does angular velocity entropy distinguish TICs from normal stem cells in practice?
- How sensitive is the metric to RNA velocity estimation noise on low-quality / multi-cell-per-spot data like Visium?

This tool is provided as-is for exploration and hypothesis generation. Interpret results in the context of your specific dataset and biological question.

## License

GNU GPLv3
