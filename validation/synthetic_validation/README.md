# Synthetic Validation

**Status:** Complete. Covered by `tests/test_velocity_synthetic.py`; figures from `make_figures.py`.

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
| Rotation      | coherent flow (±5° noise) at headings 0–180°      | ~0 at every heading | median < 0.15; rotated field gives identical scores |
| Slow noise    | fast coherent cells + slow random cells           | `min_speed` removes the slow cells | fast cells score 0, slow cells NaN |

Bifurcation now means *pairwise* angles: within-branch pairs at 0, between-branch pairs at 90°, so two of eight bins are used and the score is 1/3.

Higher dimensions (`basis='pca'`):

- 30-D aligned scores ~0; 30-D isotropic random scores ≈ the analytic `random_baseline` (0.40 for 8 bins), not 1, because random high-dimensional directions concentrate near pi/2.
- Velocities spread within a plane inside 30-D score ~1 and never exceed 1. (Under the old normalization, which divided by the random baseline, they scored ~2.5.)
- A 2D bifurcation zero-padded into 30-D gives scores identical to the 2D run.

`rotation_invariance.png` shows the fix for the old absolute-angle binning in 2D: before 0.2, the same coherent flow scored 0 at a 22.5° heading and 0.33 at 0° or 45°, because flows on a bin edge were split across two bins.

## Figures

- `synthetic_scenarios.png` — the three 2D scenarios, velocity arrows + per-cell score
- `score_distributions.png` — histograms of per-cell entropy for aligned vs bifurcation vs random
- `mixed_dataset.png` — the three-cluster mixed dataset with per-cell scores
- `rotation_invariance.png` — median score of a coherent flow vs its heading

## What this does and does not show

- Shows: the histogram + normalization math produces the expected numbers on toy input, and the score is rotation-invariant and identical across dimensions for the same configuration.
- Does not show: anything about biology, real velocity noise, or whether the metric is useful on any specific dataset. For that see the pancreas and spatial chicken heart folders.
