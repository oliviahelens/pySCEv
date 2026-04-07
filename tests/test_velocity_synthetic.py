"""
Tier 1: Synthetic validation of angular velocity entropy.

Constructs AnnData objects with known velocity fields and verifies
that score_angular_velocity_entropy produces expected results:

1. Aligned cluster: all velocity vectors point the same direction → ~0 entropy
2. Random cluster: velocity vectors uniformly random → ~1 entropy (normalized)
3. Bifurcation: two sub-clusters diverging at a branch point → intermediate entropy
4. Zero-velocity cells → NaN
5. Mixed dataset combining all of the above
"""

import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix
import pytest

from pysce._velocity import score_angular_velocity_entropy


def _make_adata(n_cells, positions, velocities):
    """Build a minimal AnnData with embedding and velocity projection."""
    X = csr_matrix((n_cells, 10))  # dummy expression matrix
    adata = AnnData(X)
    adata.obsm['X_umap'] = np.array(positions, dtype=np.float64)
    adata.obsm['velocity_umap'] = np.array(velocities, dtype=np.float64)
    return adata


# ---------------------------------------------------------------------------
# Test 1: Perfectly aligned velocities → entropy ≈ 0
# ---------------------------------------------------------------------------
def test_aligned_cluster():
    """All neighbors point the same direction — entropy should be near zero."""
    np.random.seed(42)
    n = 200
    # Tight cluster so everyone is everyone's neighbor
    positions = np.random.normal(0, 0.1, (n, 2))
    # All velocity vectors point right (+x)
    velocities = np.column_stack([np.ones(n), np.zeros(n)])

    adata = _make_adata(n, positions, velocities)
    score_angular_velocity_entropy(adata, n_neighbors=30, n_bins=8)

    scores = adata.obs['angular_velocity_entropy'].values
    valid = scores[np.isfinite(scores)]
    assert len(valid) > 0, "Expected some scored cells"
    # All vectors identical → all land in one bin → entropy = 0
    assert np.all(valid < 0.05), (
        f"Aligned cluster should have ~0 entropy, got max={valid.max():.4f}"
    )


# ---------------------------------------------------------------------------
# Test 2: Uniformly random velocities → entropy ≈ 1 (normalized)
# ---------------------------------------------------------------------------
def test_random_cluster():
    """Random directions — entropy should be near max (1.0 when normalized)."""
    np.random.seed(123)
    n = 2000  # need enough cells for histogram to fill all bins
    positions = np.random.normal(0, 0.1, (n, 2))
    # Random angles uniformly in [0, 2pi)
    theta = np.random.uniform(0, 2 * np.pi, n)
    velocities = np.column_stack([np.cos(theta), np.sin(theta)])

    adata = _make_adata(n, positions, velocities)
    score_angular_velocity_entropy(adata, n_neighbors=50, n_bins=8)

    scores = adata.obs['angular_velocity_entropy'].values
    valid = scores[np.isfinite(scores)]
    assert len(valid) > 0
    median = np.median(valid)
    assert median > 0.85, (
        f"Random cluster should have high entropy, got median={median:.4f}"
    )


# ---------------------------------------------------------------------------
# Test 3: Bifurcation — two diverging streams → intermediate entropy
# ---------------------------------------------------------------------------
def test_bifurcation():
    """
    Branch point: cells in a tight cluster with velocities pointing in
    two opposite directions (half go up-right, half go down-right).
    Entropy should be intermediate — higher than aligned, lower than random.
    """
    np.random.seed(99)
    n = 400
    positions = np.random.normal(0, 0.1, (n, 2))
    # Half point at +45°, half at -45°
    angle_a = np.pi / 4
    angle_b = -np.pi / 4
    velocities = np.zeros((n, 2))
    velocities[:n // 2] = [np.cos(angle_a), np.sin(angle_a)]
    velocities[n // 2:] = [np.cos(angle_b), np.sin(angle_b)]

    adata = _make_adata(n, positions, velocities)
    score_angular_velocity_entropy(adata, n_neighbors=50, n_bins=8)

    scores = adata.obs['angular_velocity_entropy'].values
    valid = scores[np.isfinite(scores)]
    median = np.median(valid)
    # Two bins occupied out of 8 → H = log2(2)/log2(8) = 1/3 ≈ 0.333
    # With some noise from neighbor sampling, expect roughly 0.2 - 0.5
    assert 0.1 < median < 0.6, (
        f"Bifurcation should have intermediate entropy, got median={median:.4f}"
    )


# ---------------------------------------------------------------------------
# Test 4: Zero-velocity cells → NaN
# ---------------------------------------------------------------------------
def test_zero_velocity_gives_nan():
    """Cells with zero velocity vectors should get NaN scores."""
    np.random.seed(7)
    n = 50
    positions = np.random.normal(0, 0.1, (n, 2))
    velocities = np.zeros((n, 2))  # all zero

    adata = _make_adata(n, positions, velocities)
    score_angular_velocity_entropy(adata, n_neighbors=10, n_bins=8)

    scores = adata.obs['angular_velocity_entropy'].values
    assert np.all(np.isnan(scores)), "All-zero velocities should yield all NaN"


# ---------------------------------------------------------------------------
# Test 5: Mixed dataset — aligned, random, and zero regions
# ---------------------------------------------------------------------------
def test_mixed_dataset():
    """
    Three spatial clusters: aligned (low entropy), random (high entropy),
    and zero-velocity (NaN). Verify ordering and NaN assignment.
    """
    np.random.seed(55)
    n_per = 200

    # Cluster A: aligned, centered at (-5, 0)
    pos_a = np.random.normal([-5, 0], 0.1, (n_per, 2))
    vel_a = np.tile([1.0, 0.0], (n_per, 1))

    # Cluster B: random, centered at (5, 0)
    pos_b = np.random.normal([5, 0], 0.1, (n_per, 2))
    theta = np.random.uniform(0, 2 * np.pi, n_per)
    vel_b = np.column_stack([np.cos(theta), np.sin(theta)])

    # Cluster C: zero velocity, centered at (0, 5)
    pos_c = np.random.normal([0, 5], 0.1, (n_per, 2))
    vel_c = np.zeros((n_per, 2))

    positions = np.vstack([pos_a, pos_b, pos_c])
    velocities = np.vstack([vel_a, vel_b, vel_c])

    adata = _make_adata(3 * n_per, positions, velocities)
    score_angular_velocity_entropy(adata, n_neighbors=30, n_bins=8)

    scores = adata.obs['angular_velocity_entropy'].values
    aligned_scores = scores[:n_per]
    random_scores = scores[n_per:2 * n_per]
    zero_scores = scores[2 * n_per:]

    # Aligned should be low
    aligned_valid = aligned_scores[np.isfinite(aligned_scores)]
    assert np.median(aligned_valid) < 0.1, (
        f"Aligned cluster entropy too high: {np.median(aligned_valid):.4f}"
    )

    # Random should be high
    random_valid = random_scores[np.isfinite(random_scores)]
    assert np.median(random_valid) > 0.8, (
        f"Random cluster entropy too low: {np.median(random_valid):.4f}"
    )

    # Zero should be NaN
    assert np.all(np.isnan(zero_scores)), "Zero-velocity cluster should be all NaN"

    # Ordering: aligned < random
    assert np.median(aligned_valid) < np.median(random_valid), (
        "Aligned cluster should have lower entropy than random cluster"
    )


# ---------------------------------------------------------------------------
# Test 6: .uns metadata is populated correctly
# ---------------------------------------------------------------------------
def test_metadata_stored():
    """Verify run parameters are stored in adata.uns."""
    np.random.seed(0)
    n = 100
    positions = np.random.normal(0, 0.1, (n, 2))
    velocities = np.column_stack([np.ones(n), np.zeros(n)])

    adata = _make_adata(n, positions, velocities)
    score_angular_velocity_entropy(
        adata, n_neighbors=15, n_bins=12, normalize=False, key_added='test_ent'
    )

    params = adata.uns['test_ent_params']
    assert params['basis'] == 'umap'
    assert params['n_neighbors'] == 15
    assert params['n_bins'] == 12
    assert params['normalize'] is False
    assert params['n_cells_scored'] + params['n_cells_nan'] == n


# ---------------------------------------------------------------------------
# Test 7: inplace=False returns copy without mutating original
# ---------------------------------------------------------------------------
def test_inplace_false():
    """inplace=False should return a new AnnData, leaving original untouched."""
    np.random.seed(0)
    n = 50
    positions = np.random.normal(0, 0.1, (n, 2))
    velocities = np.column_stack([np.ones(n), np.zeros(n)])

    adata = _make_adata(n, positions, velocities)
    result = score_angular_velocity_entropy(adata, inplace=False)

    assert result is not adata, "Should return a new object"
    assert 'angular_velocity_entropy' in result.obs.columns
    assert 'angular_velocity_entropy' not in adata.obs.columns


# ---------------------------------------------------------------------------
# Test 8: High-dimensional (PCA-like) embedding
# ---------------------------------------------------------------------------
def test_high_dim_pairwise():
    """
    Verify high-dim path works and gives sensible results.
    Use a 30-dim embedding with aligned vs random velocity vectors.
    """
    np.random.seed(42)
    n = 300
    ndim = 30

    # --- Aligned: all velocity vectors point along dim 0 ---
    pos_a = np.random.normal(0, 0.1, (n, ndim))
    vel_a = np.zeros((n, ndim))
    vel_a[:, 0] = 1.0

    adata_a = AnnData(csr_matrix((n, 10)))
    adata_a.obsm['X_pca'] = pos_a
    adata_a.obsm['velocity_pca'] = vel_a

    score_angular_velocity_entropy(adata_a, basis='pca', n_neighbors=30, n_bins=8)
    aligned_scores = adata_a.obs['angular_velocity_entropy'].values
    aligned_valid = aligned_scores[np.isfinite(aligned_scores)]

    # --- Random: velocity vectors in random directions ---
    pos_b = np.random.normal(0, 0.1, (n, ndim))
    vel_b = np.random.randn(n, ndim)
    vel_b = vel_b / np.linalg.norm(vel_b, axis=1, keepdims=True)

    adata_b = AnnData(csr_matrix((n, 10)))
    adata_b.obsm['X_pca'] = pos_b
    adata_b.obsm['velocity_pca'] = vel_b

    score_angular_velocity_entropy(adata_b, basis='pca', n_neighbors=30, n_bins=8)
    random_scores = adata_b.obs['angular_velocity_entropy'].values
    random_valid = random_scores[np.isfinite(random_scores)]

    assert np.median(aligned_valid) < 0.05, (
        f"Aligned high-dim entropy too high: {np.median(aligned_valid):.4f}"
    )
    # In high dimensions, pairwise angles between random unit vectors
    # concentrate around pi/2 (concentration of measure), so entropy
    # won't reach 1.0. But it should still be well above aligned (~0).
    assert np.median(random_valid) > 0.3, (
        f"Random high-dim entropy too low: {np.median(random_valid):.4f}"
    )
    assert np.median(aligned_valid) < np.median(random_valid)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
