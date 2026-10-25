# pySCEv

A fork of [pySCE](https://github.com/dchary/pysce) (Single Cell Entropy Scoring in Python) that adds an experimental RNA velocity entropy module and tests whether it provides useful biological signal — particularly in cases where expression entropy alone cannot distinguish between cell states (e.g., tumor-initiating cells vs. normal stem cells).

The original `pysce.score` expression entropy tool is unchanged. Everything under `_velocity.py` is new.

## What this fork adds

**Angular Velocity Entropy** (`pysce.score_angular_velocity_entropy`) — *Experimental.* For each cell, measures how coherent or disordered the RNA velocity vectors are in its local neighborhood. The hypothesis is that this may capture dynamic differences that expression entropy cannot: two cell populations with identical expression entropy could have very different velocity coherence if one is actively differentiating (coherent flow) while the other is dynamically disorganized.

This is speculative — see [Status](#status) and `validation/` for validation results so far.

## Installation

```bash
pip install -e .              # core
pip install -e ".[velocity]"  # + scVelo, for ensure_velocity / auto-projection
```

### Requirements

- `scanpy >= 1.8.2`, `torch >= 2.0.0`, `networkx` (expression entropy, `pysce.score`)
- `scikit-learn >= 1.0`, `anndata`, `numpy`, `scipy` (velocity entropy)
- `scvelo >= 0.2.4` (optional; only needed if pySCEv has to estimate or project velocity for you)

`pysce.score` imports torch lazily, so the velocity tools work without torch installed.

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
- `adata.obs['angular_velocity_entropy']` — per-cell score in [0, 1] (normalized), or NaN for cells slower than `min_speed`
- `adata.uns['angular_velocity_entropy_params']` — run metadata: parameters, neighbors actually used, cell counts, and `random_baseline` (the score isotropic noise would get at this dimensionality)

**Interpretation:**
- **Low entropy** (~0): neighbors' velocity vectors are coherent, pointing in similar directions (active transition / committed trajectory)
- **High entropy** (~1): neighbors' velocity vectors are scattered (disordered / quiescent / multipotent neighborhood)

### Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `basis` | `'umap'` | Embedding for neighbor lookup (`'umap'` for 2D, `'pca'` for higher-dim) |
| `n_neighbors` | `30` | Number of nearest neighbors per cell |
| `n_bins` | `8` | Angular histogram bins |
| `normalize` | `True` | Divide by `log2(n_bins)` so scores lie in [0, 1] |
| `min_speed` | `1e-10` | Cells (and neighbors) slower than this are treated as directionless |

### Normalization and the random baseline

When `normalize=True`, raw entropy is divided by `log2(n_bins)`, the largest value any neighborhood can reach in any dimension, so scores always lie in [0, 1].

What isotropic *noise* scores depends on dimensionality. Random unit vectors in high dimensions have pairwise angles concentrated near pi/2 (concentration of measure), so they fill fewer bins. That expected noise level is stored as `random_baseline` in the run metadata:

| Dimensions | Random baseline (bits, n_bins=8) | Normalized |
|---|---|---|
| 2 (UMAP) | 3.00 | 1.00 |
| 10 | 1.89 | 0.63 |
| 30 (PCA) | 1.21 | 0.40 |
| 50 | 1.05 | 0.35 |

Compare absolute scores only between runs on the same basis. (Earlier versions divided by this baseline instead; that let structured high-dimensional fields score well above 1 and was removed.)

### Velocity speed

Slow cells have poorly determined velocity directions, so their angles are dominated by estimation noise and push entropy up. On pancreas, entropy correlates with log speed (see `validation/pancreas_endocrinogenesis/`). Setting `min_speed` to a low quantile of `np.linalg.norm(adata.obsm['velocity_umap'], axis=1)` removes slow cells from every neighborhood; check that conclusions hold with it.

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
For each cell `i`, take its RNA velocity projected into embedding space (e.g., UMAP or PCA) as a vector `v_i`. Find the `k` nearest neighbors in embedding space (default `k=30`), excluding the cell itself.

**Step 2: Compute pairwise angles.**
Compute every pairwise angle among the `k` neighbor velocity vectors:

```
theta_jl = arccos( (v_j . v_l) / (|v_j| * |v_l|) ),   theta in [0, pi]
```

This gives `k*(k-1)/2` angles. They depend only on the vectors' relative orientation, so the score is the same however the embedding is rotated, and the same configuration scores identically in 2D or embedded in more dimensions.

(Versions before 0.2 binned each neighbor's *absolute* angle `atan2(v_y, v_x)` in 2D. Fixed bin edges split any flow pointing along a bin boundary — 0°, 45°, 90°, … — so an identical coherent flow could score anywhere from 0 to 0.33 depending on its heading in the UMAP.)

**Step 3: Bin and compute Shannon entropy.**
Histogram the angles into `B` equal-width bins on [0, pi] (default `B=8`). Compute the Shannon entropy of the resulting distribution:

```
p_b = count_b / sum(counts)
H = -sum_b p_b * log2(p_b)
```

If all neighbors point the same way, every angle is ~0 and `H = 0`. Two opposite streams put roughly half the pairs at 0 and half at pi, so `H ≈ 1` bit. Directions spread uniformly around a circle give uniform pairwise angles and `H = log2(B)`.

**Step 4: Normalize.**
Divide by `log2(B)`. For reference, the expected pairwise-angle density for isotropic random vectors in `d` dimensions is

```
f(theta) proportional to sin^(d-2)(theta),   theta in [0, pi]
```

which is uniform in 2D (baseline = 1.0) and concentrates near pi/2 as `d` grows (baseline ≈ 0.40 in 30D). This is reported as `random_baseline`, not used as the denominator.

**Cells slower than `min_speed`** get NaN and are left out of their neighbors' angle sets — direction is undefined when there is no motion.

## Examples

See `validation/` for figures:

- **`validation/synthetic_validation/`** — Synthetic data where the ground truth is known. Aligned velocities score ~0, bifurcation ~0.33, random ~1.0, and scores do not depend on embedding rotation.
- **`validation/pancreas_endocrinogenesis/`** — scVelo pancreas endocrinogenesis dataset (Bastidas-Ponce et al. 2019), with velocity-speed and permutation controls.
- **`validation/pancreas_endocrinogenesis_enrichment/`** — gene programs of low- vs high-entropy cells on pancreas.
- **`validation/spatial_chicken_heart/`** — SIRV-imputed Visium (Mantri et al. 2021, day 14). Spatial vs UMAP neighbors. *Numbers predate the 0.2 scorer and need re-running.*
- **`validation/parameterization/`** — k neighbors × bin count sweep on pancreas. Median Spearman rho = 0.95 across settings; the two lowest cell types are stable, the top of the ranking is not.

## Status

**Angular Velocity Entropy** is speculative and under active validation. The core question is whether velocity-based entropy provides biological insight that expression entropy and existing tools (e.g., scVelo's velocity confidence) do not.

**Notes so far** (pySCEv 0.2; see "Changes in 0.2" below for what moved):
- **Synthetic:** aligned (0), bifurcation (0.33), random (1.0); identical scores under rotation and across 2D / zero-padded 30D.
- **Pancreas:** the endocrine-progenitor wave (Ngn3 high EP, Pre-endocrine) scores far lower than every other population. This survives removing the slowest 25% of cells, regressing out velocity speed, and a shuffled-velocity permutation null (Kruskal-Wallis H 1811 vs null 129 ± 65, z = 26).
- **Velocity speed is a major confound.** Entropy correlates with log speed (Spearman -0.49). Once speed is accounted for, Ductal, Ngn3 low EP and the terminal endocrine types are roughly indistinguishable; "terminal cells score higher" is mostly "slow cells score higher". Use `min_speed` and check.
- **Relation to existing metrics:** anti-correlates with scVelo confidence (r = -0.47), weakly with expression entropy (r = 0.15). It is close to the simpler mean angular deviation in rank terms (Spearman 0.94, Pearson 0.78), so the added value of the entropy form is small on this dataset.
- **Gene programs:** low entropy tracks the endocrine-progenitor TF program (Neurog3, Pax4, Nkx6-1; module rho -0.56). The earlier "high entropy = mature/secretory within each cell type" result did not reproduce: within Alpha and Beta, hormone genes are higher in the *low*-entropy cells.
- **Kernel smoothing:** about half of the spatial autocorrelation on the pancreas UMAP (Moran's I null 0.45 of observed 0.93) comes from overlapping k-NN neighborhoods alone, as the chicken heart permutation also found. Smooth patches are not evidence by themselves.
- **Parameterization is stable:** median Spearman rho = 0.95 across a 5×5 sweep of k and bin count; Ngn3 high EP and Pre-endocrine are the two lowest in all 25 settings.
- **Spatial chicken heart:** the old results (spatial vs UMAP neighbors r = 0.39; low entropy in the left ventricular wall, high in valves and atria; permutation z = 5.7 for region ranking) were computed with the pre-0.2 scorer and are pending a re-run.

**Changes in 0.2:**
- 2D scores used fixed absolute-angle bins, so a coherent flow pointing along a bin edge scored ~0.33 instead of ~0 (heading-dependent). Scoring now uses pairwise angles in every dimension, which is rotation-invariant.
- High-dimensional scores were divided by the random-vector baseline and could exceed 1 (~2.5 for a field spread in a plane). Scores are now divided by log2(n_bins) and always lie in [0, 1]; the baseline is reported as metadata.
- Each cell's nearest neighbor was dropped by an off-by-one in the k-NN lookup.
- `pysce.score` crashed on NumPy 2 before scoring anything.

**Open questions:**
- Does angular velocity entropy distinguish TICs from normal stem cells in practice?
- How sensitive is the metric to RNA velocity estimation noise on low-quality / multi-cell-per-spot data like Visium?

This tool is provided as-is for exploration and hypothesis generation. Interpret results in the context of your specific dataset and biological question.

## License

GNU GPLv3
