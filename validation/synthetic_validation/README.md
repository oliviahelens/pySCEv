# Synthetic Validation

**Status:** Complete. Covered by `tests/test_velocity_synthetic.py`.

## Why

Check that `score_angular_velocity_entropy` returns the numbers it should when the velocity field is constructed by hand. This is a math/implementation check, not a biology claim.

## Scenarios

All cells placed in a tight 2D cluster so every cell sees the same neighborhood. `n_neighbors=30-50`, `n_bins=8`, normalized.

| scenario      | velocity field                                    | expected       | measured assertion                         |
|---------------|---------------------------------------------------|----------------|--------------------------------------------|
| Aligned       | every vector `(1, 0)`                             | 0              | median < 0.05                              |
| Bifurcation   | half at +45 deg, half at -45 deg                  | log2(2)/log2(8) = 1/3 | 0.1 < median < 0.6                  |
| Random        | directions uniform on the circle                  | 1.0            | median > 0.85                              |
| Zero velocity | all vectors `(0, 0)`                              | NaN            | all NaN                                    |
| Mixed         | three spatial clusters: aligned / random / zero   | ordering       | aligned < random, zero-cluster all NaN     |

Also tested in 30-D (`basis='pca'`): aligned median < 0.05, random median > 0.85. This exercises the pairwise-angle path and the dimension-aware normalization — random unit vectors in 30 D concentrate near pi/2, so without dimension-aware normalization they would score well below 1.0.

## Figures

- `synthetic_scenarios.png` — the three 2D scenarios, velocity arrows + per-cell score
- `score_distributions.png` — histograms of per-cell entropy for aligned vs bifurcation vs random
- `mixed_dataset.png` — the three-cluster mixed dataset with per-cell scores

## What this does and does not show

- Shows: the histogram + normalization math produces the expected numbers on toy input, and the high-dim path is consistent with the 2D path.
- Does not show: anything about biology, real velocity noise, or whether the metric is useful on any specific dataset. For that see the pancreas and spatial chicken heart folders.
