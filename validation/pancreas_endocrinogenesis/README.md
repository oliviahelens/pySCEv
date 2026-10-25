# Pancreas Endocrinogenesis

**Status:** Re-run on the pySCEv 0.2 scorer (pairwise angles, rotation-invariant), with speed and permutation controls added. UMAP basis, UMAP neighbors. All numbers below are in `summary.json`.

## Dataset

- scVelo built-in `pancreas` dataset (Bastidas-Ponce et al. 2019, Development)
- 3,696 cells, spliced/unspliced available
- Cell types: Ductal, Ngn3 low EP, Ngn3 high EP, Pre-endocrine, Alpha, Beta, Delta, Epsilon

## Why

First test on real scRNA-seq with a well-characterized differentiation hierarchy. The purpose is to see whether angular velocity entropy picks up any structure at all, and whether that structure is (a) interpretable, (b) non-redundant with metrics that already exist, and (c) not just a side effect of velocity speed or of the scoring kernel.

## How

```
python validation/pancreas_endocrinogenesis/run_analysis.py --permute 100
```

`validation/_common/score_pipeline.py`: scVelo preprocessing (2000 HVGs, `moments` with 30 PCs / 30 neighbors, deterministic velocity, `velocity_graph`, `velocity_embedding`), then `pysce.score_angular_velocity_entropy(adata, basis='umap', n_neighbors=30, n_bins=8)`. Expression entropy (`pysce.score`) runs on spliced+unspliced counts with mouse symbols uppercased to match the human PPI (10,063 shared genes); it takes ~80 min on CPU and is cached in `cache/`.

## Results (k=30, bins=8, UMAP basis)

### Cell-type ranking

| Cell type | Median entropy | With slowest 25% removed | Speed-residual median |
|---|---:|---:|---:|
| Ngn3 high EP | **0.00** | 0.00 | **-0.35** |
| Pre-endocrine | **0.12** | 0.07 | **-0.23** |
| Beta | 0.36 | 0.20 | +0.10 |
| Alpha | 0.60 | 0.46 | +0.20 |
| Ngn3 low EP | 0.66 | 0.31 | +0.13 |
| Ductal | 0.77 | 0.35 | +0.13 |
| Epsilon | 0.79 | 0.43 | +0.13 |
| Delta | 0.87 | 0.40 | +0.18 |

Lowest: **Ngn3 high EP and Pre-endocrine**, the populations riding the endocrine specification wave. 65% of Ngn3 high EP cells score exactly 0 (all neighbor velocity pairs within 22.5°). This is the one result that is robust to every control below.

The top of the ranking changed from the old scorer: Ductal and Epsilon now sit with Delta at the top, and Alpha is mid-pack. Among the remaining six types the order is not stable (see the speed control and `../parameterization/`).

### Relationship to other metrics

| Comparison | Pearson r | Spearman rho |
|---|---:|---:|
| scVelo velocity confidence | -0.47 | -0.51 |
| Mean angular deviation (first moment) | 0.78 | 0.94 |
| log10 velocity speed | -0.52 | -0.49 |
| Expression entropy (pySCE) | 0.15 | 0.18 |

- **Expression entropy:** weakly related (was reported as r ~ 0.02 with the old scorer). Still not a stand-in for potency.
- **Velocity confidence:** about the same as before (-0.47 vs -0.48).
- **Mean angular deviation:** Spearman 0.94. In rank terms the entropy is now nearly the same as the simpler mean-deviation score, so the earlier "~20% unique variance from multimodal neighborhoods" claim is weaker than it looked; most of that difference was binning noise.

### Control 1: velocity speed (`speed_control.png`)

Entropy falls with speed (rho -0.49): slow cells have poorly determined directions, so they score high regardless of biology. Two checks:

- **Remove the slowest 25%** (`min_speed` = 0.0122, the 25th percentile): Kruskal-Wallis H across cell types drops from 1811 to 1246. Ngn3 high EP and Pre-endocrine stay lowest; the other six compress into 0.20–0.46 and reorder (Alpha becomes highest, Ductal drops to mid). Rank correlation of cell-type medians with the raw ranking: 0.74.
- **Regress entropy on log speed** and look at residuals: H = 1618. Ngn3 high EP (-0.35) and Pre-endocrine (-0.23) are still clearly below zero; the other six all sit at +0.10 to +0.20 with overlapping boxes. Rank correlation with the raw ranking: 0.69.

**Read:** "progenitor wave = coherent flow" is not explained by speed. The separation *among* Ductal, Ngn3 low EP and the terminal endocrine types is mostly speed.

### Control 2: permutation (`permutation_null.png`)

Shuffle the UMAP velocity vectors across cells (N=100) and rescore:

| Statistic | Observed | Null mean ± sd | z |
|---|---:|---:|---:|
| Kruskal-Wallis H (cell type) | 1811 | 129 ± 65 | 25.9 |
| Moran's I (UMAP 30-NN graph) | 0.93 | 0.45 ± 0.04 | 12.2 |

Shuffled fields score ~0.99 in every cell type. The cell-type separation is far beyond what the kernel produces on its own. About half of the spatial smoothness (Moran's I null 0.45 of observed 0.93) comes from overlapping k-NN neighborhoods alone, matching what the chicken heart permutation found. Smooth patches on the UMAP are not evidence by themselves.

## Figures

- `umap_panels.png` — UMAP colored by cell type + angular velocity entropy
- `velocity_streamlines.png` — scVelo streamlines over the entropy field
- `entropy_by_celltype.png` — box plot of entropy per cell type
- `entropy_vs_confidence.png` — scatter vs scVelo velocity confidence
- `entropy_vs_angular_deviation.png` — entropy vs mean angular deviation
- `expression_vs_velocity_entropy.png` — UMAPs comparing expression entropy vs velocity entropy
- `expression_vs_velocity_scatter.png` — per-cell scatter of the two entropies
- `velocity_vs_entropy.png` — entropy vs velocity magnitude (from `velocity_vs_entropy.py`)
- `speed_control.png` — entropy vs speed, cell-type medians with slow cells removed, speed-residual box plot
- `permutation_null.png` — observed Kruskal-Wallis H and Moran's I against the shuffled-velocity null
- `alpha_phase_portraits.png` — phase portraits for Alpha-driver genes (velocity sanity check; independent of the entropy metric, not regenerated)

## Caveats

- **One dataset, one velocity model.** Deterministic scVelo mode; parameter stability is covered in `validation/parameterization/`.
- **UMAP-projected velocity.** `velocity_umap` is a kernel-smoothed projection of high-dimensional velocity, so neighborhood coherence is partly inherited from that smoothing. The permutation null shuffles projected vectors and so does not capture this. A PCA-basis run would help.
- **Ties at 0.** 16% of all cells score exactly 0, which matters for any quantile-based downstream analysis.
- The correlations above are from a single run; not bootstrapped.

## What this does and does not show

- **Shows:** the endocrine-progenitor wave (Ngn3 high EP → Pre-endocrine) has far more coherent local velocity than everything else. This holds after removing slow cells, after regressing out speed, and against a shuffled-velocity null.
- **Does not show:** a meaningful ordering among Ductal, Ngn3 low EP and the terminal endocrine types. Those differences are largely speed. Also does not show that the metric adds much over mean angular deviation (Spearman 0.94), or that it generalizes beyond this dataset.
