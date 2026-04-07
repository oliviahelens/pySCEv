####################################################################################################
# # Copyright (C) 2024-Present - Daniel Charytonowicz - All Rights Reserved
# # Contact: daniel.charytonowicz@icahn.mssm.edu
# ###################################################################################################

from typing import Optional
from anndata import AnnData
import numpy as np


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
    key_added: str = 'angular_velocity_entropy',
    inplace: bool = True,
) -> Optional[AnnData]:
    """\
    Local Angular Velocity Entropy (Metric 1)

    For each cell, examines the RNA velocity vectors of its k nearest
    neighbors in embedding space and computes Shannon entropy of how
    spread out those velocity directions are among the neighbors.

    High entropy = neighbors moving in scattered directions
    (disordered / multipotent neighborhood).
    Low entropy = neighbors moving coherently in the same direction
    (committed trajectory).

    For 2D embeddings, each neighbor's velocity direction is converted
    to an absolute angle via atan2 and binned in [0, 2pi). For higher-
    dimensional embeddings, all pairwise angles among neighbor velocity
    vectors are computed via cosine similarity and binned in [0, pi].

    Cells with near-zero velocity (magnitude < 1e-10) are assigned
    NaN. Summary statistics (number of scored cells, NaN count,
    parameters used) are stored in adata.uns[key_added + '_params'].

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
        If True, normalize entropy to [0, 1] by dividing by
        log2(n_bins) (the maximum possible entropy).
    key_added
        Key in adata.obs where scores are stored.
    inplace
        Whether to modify adata in-place.

    Returns
    -------
    AnnData or None if inplace.
    """
    from sklearn.neighbors import NearestNeighbors
    from tqdm.auto import tqdm

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

    V = np.asarray(adata.obsm[velocity_key])
    X_emb = np.asarray(adata.obsm[embedding_key])
    n_cells = adata.n_obs
    ndim = V.shape[1]

    # Identify cells with near-zero velocity (undefined direction)
    vel_norms = np.linalg.norm(V, axis=1)
    zero_vel_mask = vel_norms < 1e-10

    # Build k-NN index in embedding space
    k = min(n_neighbors, n_cells - 1)
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm='auto')
    nn.fit(X_emb)
    indices = nn.kneighbors(return_distance=False)[:, 1:]  # exclude self

    entropy_scores = np.full(n_cells, np.nan)

    if ndim == 2:
        # 2D: absolute angle of each velocity vector, binned in [0, 2pi).
        # Entropy measures how spread out neighbor directions are.
        angles = np.arctan2(V[:, 1], V[:, 0]) % (2 * np.pi)
        bin_edges = np.linspace(0, 2 * np.pi, n_bins + 1)

        for i in tqdm(range(n_cells), desc="Scoring angular velocity entropy", unit="cell"):
            if zero_vel_mask[i]:
                continue
            nbr_idx = indices[i][~zero_vel_mask[indices[i]]]
            if len(nbr_idx) == 0:
                continue
            counts, _ = np.histogram(angles[nbr_idx], bins=bin_edges)
            p = counts[counts > 0] / counts.sum()
            entropy_scores[i] = -np.sum(p * np.log2(p))
    else:
        # Higher-dim: pairwise angles among neighbor velocity vectors
        # via cosine similarity, binned in [0, pi]. Measures the same
        # thing as the 2D path — directional spread — but generalized
        # to arbitrary dimensions.
        V_unit = np.zeros_like(V)
        valid = ~zero_vel_mask
        V_unit[valid] = V[valid] / vel_norms[valid, np.newaxis]
        bin_edges = np.linspace(0, np.pi, n_bins + 1)

        for i in tqdm(range(n_cells), desc="Scoring angular velocity entropy", unit="cell"):
            if zero_vel_mask[i]:
                continue
            nbr_idx = indices[i][~zero_vel_mask[indices[i]]]
            if len(nbr_idx) < 2:
                continue
            nbr_vecs = V_unit[nbr_idx]
            # Pairwise cosine similarities (upper triangle)
            cos_matrix = nbr_vecs @ nbr_vecs.T
            triu_idx = np.triu_indices(len(nbr_idx), k=1)
            cos_pairs = np.clip(cos_matrix[triu_idx], -1.0, 1.0)
            theta = np.arccos(cos_pairs)
            counts, _ = np.histogram(theta, bins=bin_edges)
            p = counts[counts > 0] / counts.sum()
            entropy_scores[i] = -np.sum(p * np.log2(p))

    # Normalize to [0, 1] if requested
    if normalize:
        max_entropy = np.log2(n_bins)
        if max_entropy > 0:
            entropy_scores = entropy_scores / max_entropy

    adata.obs[key_added] = entropy_scores

    # Store run metadata so users can inspect parameters and coverage
    n_scored = int(np.isfinite(entropy_scores).sum())
    n_nan = int(np.isnan(entropy_scores).sum())
    adata.uns[key_added + '_params'] = {
        'basis': basis,
        'n_neighbors': n_neighbors,
        'n_bins': n_bins,
        'normalize': normalize,
        'n_cells_scored': n_scored,
        'n_cells_nan': n_nan,
    }

    if not inplace:
        return adata
