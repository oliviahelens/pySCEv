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
    assert params['n_neighbors_used'] == 15
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
    # Random directions in 30D concentrate near pi/2, so they should land
    # near the analytic random baseline rather than at the [0, 1] ceiling.
    baseline = adata_b.uns['angular_velocity_entropy_params']['random_baseline']
    assert abs(np.median(random_valid) - baseline) < 0.05, (
        f"Random high-dim entropy {np.median(random_valid):.4f} "
        f"far from baseline {baseline:.4f}"
    )
    assert np.median(aligned_valid) < np.median(random_valid)


# ---------------------------------------------------------------------------
# Test 9: Score does not depend on how the embedding is rotated
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('heading_deg', [0, 10, 22.5, 45, 180])
def test_rotation_invariance(heading_deg):
    """
    A coherent flow with small angular noise must score the same whatever
    direction it points. Fixed absolute-angle bins used to split flows
    that sat on a bin edge (0, 45, 90 ... degrees) and score them ~0.33.
    """
    rng = np.random.default_rng(1)
    n = 400
    positions = rng.normal(0, 0.1, (n, 2))
    theta = np.deg2rad(heading_deg) + rng.normal(0, np.deg2rad(5), n)
    velocities = np.column_stack([np.cos(theta), np.sin(theta)])

    adata = _make_adata(n, positions, velocities)
    score_angular_velocity_entropy(adata, n_neighbors=30, n_bins=8)
    median = np.median(adata.obs['angular_velocity_entropy'])
    assert median < 0.15, f"heading {heading_deg}: median={median:.3f}"


def test_rotation_gives_identical_scores():
    """Rotating every velocity vector by the same angle changes nothing."""
    rng = np.random.default_rng(2)
    n = 300
    positions = rng.normal(0, 1, (n, 2))
    velocities = rng.normal(0, 1, (n, 2)) + [2, 0]
    a = np.deg2rad(37)
    R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])

    ad1 = _make_adata(n, positions, velocities)
    ad2 = _make_adata(n, positions, velocities @ R.T)
    score_angular_velocity_entropy(ad1)
    score_angular_velocity_entropy(ad2)
    np.testing.assert_allclose(
        ad1.obs['angular_velocity_entropy'], ad2.obs['angular_velocity_entropy'],
        atol=1e-12,
    )


# ---------------------------------------------------------------------------
# Test 10: Normalized scores never exceed 1, in any dimension
# ---------------------------------------------------------------------------
def test_high_dim_bounded():
    """
    Velocities spread uniformly within a 2D plane of a 30D space give
    pairwise angles uniform on [0, pi] -- more spread than isotropic noise.
    These used to score ~2.5 when normalized by the random baseline.
    """
    rng = np.random.default_rng(3)
    n, ndim = 400, 30
    theta = rng.uniform(0, 2 * np.pi, n)
    V = np.zeros((n, ndim))
    V[:, 0], V[:, 1] = np.cos(theta), np.sin(theta)

    adata = AnnData(csr_matrix((n, 10)))
    adata.obsm['X_pca'] = rng.normal(0, 0.1, (n, ndim))
    adata.obsm['velocity_pca'] = V
    score_angular_velocity_entropy(adata, basis='pca')
    scores = adata.obs['angular_velocity_entropy'].values
    assert np.nanmax(scores) <= 1.0 + 1e-12
    assert np.median(scores) > 0.85


# ---------------------------------------------------------------------------
# Test 11: Same configuration scores the same in 2D and embedded in 30D
# ---------------------------------------------------------------------------
def test_bifurcation_consistent_across_dims():
    """Pairwise angles are dimension-free, so a 2D bifurcation embedded in
    30D (zero padding) must give identical scores."""
    rng = np.random.default_rng(4)
    n = 300
    sign = np.where(rng.random(n) < 0.5, 1.0, -1.0)
    V2 = np.column_stack([sign, np.zeros(n)]) + rng.normal(0, 0.05, (n, 2))
    pos2 = rng.normal(0, 0.1, (n, 2))

    ad2 = _make_adata(n, pos2, V2)
    ad30 = AnnData(csr_matrix((n, 10)))
    ad30.obsm['X_pca'] = np.pad(pos2, ((0, 0), (0, 28)))
    ad30.obsm['velocity_pca'] = np.pad(V2, ((0, 0), (0, 28)))
    score_angular_velocity_entropy(ad2)
    score_angular_velocity_entropy(ad30, basis='pca')
    np.testing.assert_allclose(
        ad2.obs['angular_velocity_entropy'], ad30.obs['angular_velocity_entropy'],
        atol=1e-12,
    )
    # Opposite streams: ~half the pairs at 0, half at pi -> ~1 bit of 3
    assert 0.25 < np.median(ad2.obs['angular_velocity_entropy']) < 0.4


# ---------------------------------------------------------------------------
# Test 12: min_speed drops slow cells and slow neighbors
# ---------------------------------------------------------------------------
def test_min_speed_excludes_slow_noise():
    """
    Coherent fast cells mixed with slow cells pointing randomly. With the
    default threshold the slow cells' noise inflates entropy; raising
    min_speed above their speed removes them from every neighborhood.
    """
    rng = np.random.default_rng(5)
    n = 400
    positions = rng.normal(0, 0.1, (n, 2))
    fast = np.arange(n) < n // 2
    theta = np.where(fast, 0.3, rng.uniform(0, 2 * np.pi, n))
    speed = np.where(fast, 1.0, 1e-3)
    velocities = speed[:, None] * np.column_stack([np.cos(theta), np.sin(theta)])

    adata = _make_adata(n, positions, velocities)
    score_angular_velocity_entropy(adata, key_added='all')
    score_angular_velocity_entropy(adata, min_speed=1e-2, key_added='fast')

    assert np.median(adata.obs['all'][fast]) > 0.5
    assert np.all(adata.obs['fast'][fast] < 1e-12)
    assert np.all(np.isnan(adata.obs['fast'][~fast]))
    assert adata.uns['fast_params']['min_speed'] == 1e-2


# ---------------------------------------------------------------------------
# Test 13: Vectorized implementation matches a direct per-cell reference
# ---------------------------------------------------------------------------
def test_matches_reference_loop():
    """Cross-check chunked vectorized scoring against a plain loop."""
    from pysce._velocity import _pairwise_angle_entropy
    rng = np.random.default_rng(6)
    n, d, k, B = 120, 5, 12, 8
    V = rng.normal(size=(n, d))
    V[::7] = 0  # some invalid cells
    norms = np.linalg.norm(V, axis=1)
    valid = norms > 1e-10
    U = np.zeros_like(V)
    U[valid] = V[valid] / norms[valid, None]
    idx = np.array([rng.choice(np.delete(np.arange(n), i), k, replace=False)
                    for i in range(n)])

    got = _pairwise_angle_entropy(U, valid, idx, B, chunk_size=17)

    expected = np.full(n, np.nan)
    for i in range(n):
        nb = idx[i][valid[idx[i]]]
        if len(nb) < 2:
            continue
        c = np.clip(U[nb] @ U[nb].T, -1, 1)[np.triu_indices(len(nb), 1)]
        counts, _ = np.histogram(np.arccos(c), bins=np.linspace(0, np.pi, B + 1))
        p = counts[counts > 0] / counts.sum()
        expected[i] = -np.sum(p * np.log2(p))
    np.testing.assert_allclose(got, expected, atol=1e-10)


def test_import_without_torch():
    """The velocity tools must not require torch."""
    import subprocess, sys
    code = (
        "import sys; sys.modules['torch'] = None; "
        "import pysce; pysce.score_angular_velocity_entropy"
    )
    subprocess.run([sys.executable, '-c', code], check=True)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
