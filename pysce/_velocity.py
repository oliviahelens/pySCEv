####################################################################################################
# # Copyright (C) 2024-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from typing import Optional
from anndata import AnnData
import numpy as np


def _random_baseline_entropy(ndim, n_bins):
    """\
    Expected Shannon entropy (in bits) of the binned pairwise-angle
    distribution for isotropic random unit vectors in ``ndim`` dimensions,
    in the limit of many neighbors.

    Pairwise angles between random unit vectors follow
    f(theta) ∝ sin^(d-2)(theta) on [0, pi]. For d = 2 this is uniform, so
    the baseline equals the true maximum log2(n_bins). For larger d the
    density concentrates near pi/2 (concentration of measure), so random
    vectors fill fewer bins and the baseline drops well below the maximum.

    This is a *reference value*, not an upper bound: structured
    configurations (e.g. velocities spread within a low-dimensional
    subspace) can exceed it. It is reported in the run metadata so scores
    can be read against what pure noise would give at that dimensionality.

    Reference values (n_bins=8):

        ndim    baseline (bits)    as fraction of log2(8)
        ----    ---------------    ----------------------
           2              3.00        1.00
          10              1.89        0.63
          30              1.21        0.40
          50              1.05        0.35

    Params
    -------
    ndim
        Dimensionality of the velocity vectors.
    n_bins
        Number of angular histogram bins.

    Returns
    -------
    baseline : float
        Expected Shannon entropy (bits) for random directions.
    """
    if ndim == 2:
        return float(np.log2(n_bins))

    from scipy.integrate import quad

    bin_edges = np.linspace(0, np.pi, n_bins + 1)
    bin_masses = np.array([
        quad(lambda t: np.sin(t) ** (ndim - 2), bin_edges[j], bin_edges[j + 1])[0]
        for j in range(n_bins)
    ])
    p = bin_masses / bin_masses.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def _pairwise_angle_entropy(V_unit, valid, indices, n_bins, chunk_size=2048):
    """\
    Shannon entropy (bits) of the binned pairwise angles among each cell's
    neighbor velocity directions.

    Pairwise angles depend only on the vectors' relative orientation, so
    the result is invariant to rotations of the embedding. Neighbors whose
    velocity is flagged invalid are dropped from every pair they take part
    in. Cells with fewer than two valid neighbors get NaN.
    """
    n_cells, k = indices.shape
    triu_j, triu_l = np.triu_indices(k, k=1)
    out = np.full(n_cells, np.nan)

    for start in range(0, n_cells, chunk_size):
        idx = indices[start:start + chunk_size]
        nbr = V_unit[idx]                                   # (c, k, d)
        nbr_valid = valid[idx]                              # (c, k)

        cos = np.einsum('ckd,cld->ckl', nbr, nbr)[:, triu_j, triu_l]
        theta = np.arccos(np.clip(cos, -1.0, 1.0))          # (c, P)
        bins = np.minimum((theta / np.pi * n_bins).astype(np.intp), n_bins - 1)
        pair_ok = nbr_valid[:, triu_j] & nbr_valid[:, triu_l]

        counts = np.zeros((len(idx), n_bins))
        rows = np.broadcast_to(np.arange(len(idx))[:, None], bins.shape)
        np.add.at(counts, (rows[pair_ok], bins[pair_ok]), 1)

        total = counts.sum(1, keepdims=True)
        with np.errstate(divide='ignore', invalid='ignore'):
            p = counts / total
            h = -np.sum(np.where(p > 0, p * np.log2(p), 0.0), axis=1)
        h[total[:, 0] == 0] = np.nan
        out[start:start + len(idx)] = h

    return out


def ensure_velocity(
    adata: AnnData,
    mode: str = 'stochastic',
    basis: str = 'umap',
) -> None:
    """\
    Ensure Velocity

    Ensures that adata contains all velocity components needed for
    downstream entropy scoring: velocity layer, velocity graph, and
    velocity embedding projection.

    Runs each step only if the corresponding output is missing.

    .. warning::
        When velocity needs to be estimated from scratch, this calls
        ``scv.pp.filter_and_normalize`` which **subsets genes**
        (default ``n_top_genes=2000``) **and normalizes counts
        in-place**. If you need to preserve the original gene set
        or have custom preprocessing, run scVelo yourself on a copy
        and pass the precomputed velocity layers to
        ``score_angular_velocity_entropy`` directly.

    Params
    -------
    adata
        Annotated data matrix with spliced/unspliced layers
        (required only if velocity layer is missing).
    mode
        scVelo velocity mode: 'stochastic', 'deterministic', or 'dynamical'.
    basis
        Embedding basis for velocity projection (e.g. 'umap', 'pca').
    """
    import scvelo as scv
    import warnings

    if 'velocity' not in adata.layers:

        warnings.warn(
            "No precomputed velocity found. Running scVelo pipeline which "
            "will filter genes (n_top_genes=2000) and normalize adata "
            "in-place. To avoid this, precompute velocity on a copy and "
            "pass it directly.",
            UserWarning,
            stacklevel=2,
        )

        if 'spliced' not in adata.layers or 'unspliced' not in adata.layers:
            raise ValueError(
                "AnnData must contain 'spliced' and 'unspliced' layers "
                "for velocity estimation. Run velocyto or kb-python "
                "to generate these counts."
            )

        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

        if mode == 'dynamical':
            scv.tl.recover_dynamics(adata)

        scv.tl.velocity(adata, mode=mode)

    if 'velocity_graph' not in adata.uns:
        scv.tl.velocity_graph(adata)

    velocity_key = f'velocity_{basis}'
    if velocity_key not in adata.obsm:
        if f'X_{basis}' not in adata.obsm:
            raise ValueError(
                f"Embedding 'X_{basis}' not found in adata.obsm. "
                f"Cannot project velocity onto '{basis}'. Compute "
                f"the embedding first (e.g. sc.tl.umap)."
            )
        scv.tl.velocity_embedding(adata, basis=basis)


def score_angular_velocity_entropy(
    adata: AnnData,
    basis: str = 'umap',
    n_neighbors: int = 30,
    n_bins: int = 8,
    normalize: bool = True,
    min_speed: float = 1e-10,
    key_added: str = 'angular_velocity_entropy',
    inplace: bool = True,
) -> Optional[AnnData]:
    """\
    Local Angular Velocity Entropy

    For each cell, takes the RNA velocity vectors of its k nearest
    neighbors in embedding space, computes every pairwise angle between
    them, and returns the Shannon entropy of the binned angles.

    Low entropy = neighbors moving coherently in the same direction.
    High entropy = neighbors moving in scattered directions.

    Pairwise angles are used for every dimensionality (2D included), so
    the score does not depend on how the embedding happens to be rotated.
    Angles lie in [0, pi] and are split into ``n_bins`` equal-width bins.

    Cells whose own velocity is slower than ``min_speed`` are assigned NaN,
    and such neighbors are excluded from the pairwise angles. Slow cells
    have poorly determined directions, so their angles are dominated by
    noise; raise ``min_speed`` (e.g. to a low quantile of the speed
    distribution) to keep that noise out of the score. Run metadata is
    stored in adata.uns[key_added + '_params'].

    Params
    -------
    adata
        Annotated data matrix. Must have velocity projection in
        adata.obsm[f'velocity_{basis}']. If missing, attempts
        to compute it via scVelo.
    basis
        Embedding for neighbor lookup and velocity projection.
        'umap' for 2D, 'pca' for higher-dimensional.
    n_neighbors
        Number of nearest neighbors per cell.
    n_bins
        Number of angular bins for the histogram.
    normalize
        If True, divide by log2(n_bins), the largest entropy any
        neighborhood can have, so scores lie in [0, 1]. The expected
        score for isotropic random directions depends on dimensionality
        (1.0 in 2D, about 0.4 in 30D for 8 bins) and is stored as
        'random_baseline' in the run metadata. Compare absolute scores
        only between runs on the same basis.
    min_speed
        Velocity vectors with L2 norm below this value are treated as
        having no direction.
    key_added
        Key in adata.obs where scores are stored.
    inplace
        Whether to modify adata in-place.

    Returns
    -------
    AnnData or None if inplace.
    """
    from sklearn.neighbors import NearestNeighbors

    velocity_key = f'velocity_{basis}'
    embedding_key = f'X_{basis}'

    # Copy upfront if not inplace, so we never mutate the caller's object
    if not inplace:
        adata = adata.copy()

    # Validate embedding exists
    if embedding_key not in adata.obsm:
        raise ValueError(
            f"Embedding '{embedding_key}' not found in adata.obsm. "
            f"Run the corresponding embedding method first."
        )

    # Ensure velocity projection exists
    if velocity_key not in adata.obsm:
        if 'velocity' not in adata.layers:
            raise ValueError(
                "No velocity data found. Either precompute velocity "
                "(e.g. via scvelo or pysce.ensure_velocity) or provide "
                f"adata.obsm['{velocity_key}'] directly."
            )
        try:
            import scvelo as scv
        except ImportError:
            raise ImportError(
                f"Velocity embedding '{velocity_key}' not found and scvelo "
                "is not installed. Either precompute velocity embeddings or "
                "install scvelo: pip install scvelo"
            )
        if 'velocity_graph' not in adata.uns:
            scv.tl.velocity_graph(adata)
        scv.tl.velocity_embedding(adata, basis=basis)

    V = np.asarray(adata.obsm[velocity_key], dtype=np.float64)
    X_emb = np.asarray(adata.obsm[embedding_key])
    n_cells = adata.n_obs
    ndim = V.shape[1]

    if n_cells < 3:
        raise ValueError("Need at least 3 cells to score neighborhoods.")

    # Cells with near-zero (or non-finite) velocity have no usable direction
    vel_norms = np.linalg.norm(V, axis=1)
    valid = np.isfinite(vel_norms) & (vel_norms >= min_speed)
    V_unit = np.zeros_like(V)
    V_unit[valid] = V[valid] / vel_norms[valid, np.newaxis]

    # Build k-NN index in embedding space
    k = min(n_neighbors, n_cells - 1)
    nn = NearestNeighbors(n_neighbors=k, algorithm='auto')
    nn.fit(X_emb)
    # With no query passed, sklearn excludes each cell from its own list
    indices = nn.kneighbors(return_distance=False)

    entropy_scores = _pairwise_angle_entropy(V_unit, valid, indices, n_bins)
    entropy_scores[~valid] = np.nan

    max_entropy = float(np.log2(n_bins))
    baseline = _random_baseline_entropy(ndim, n_bins)
    if normalize:
        entropy_scores = entropy_scores / max_entropy

    adata.obs[key_added] = entropy_scores

    # Store run metadata so users can inspect parameters and coverage
    adata.uns[key_added + '_params'] = {
        'basis': basis,
        'embedding_ndim': ndim,
        'n_neighbors': n_neighbors,
        'n_neighbors_used': k,
        'n_bins': n_bins,
        'angles': 'pairwise',
        'normalize': normalize,
        'min_speed': float(min_speed),
        'max_entropy_bits': max_entropy,
        'random_baseline_bits': baseline,
        'random_baseline': baseline / max_entropy if normalize else baseline,
        'n_cells_scored': int(np.isfinite(entropy_scores).sum()),
        'n_cells_nan': int(np.isnan(entropy_scores).sum()),
    }

    if not inplace:
        return adata
