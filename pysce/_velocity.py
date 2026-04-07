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
) -> None:
    """\
    Ensure Velocity

    Runs scVelo velocity estimation if 'velocity' layer is missing.
    Expects adata to contain 'spliced' and 'unspliced' layers.
    Modifies adata in-place, adding velocity layers, velocity graph,
    and velocity embedding.

    Params
    -------
    adata
        Annotated data matrix with spliced/unspliced layers.
    mode
        scVelo velocity mode: 'stochastic', 'deterministic', or 'dynamical'.
    """
    import scvelo as scv

    if 'velocity' not in adata.layers:

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

    # Ensure velocity graph exists (needed for velocity_embedding)
    if 'velocity_graph' not in adata.uns:
        scv.tl.velocity_graph(adata)


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
    neighbors in embedding space and computes Shannon entropy of the
    angular distribution.

    High entropy = disordered, multipotent neighborhood dynamics.
    Low entropy = coherent directional flow (committed trajectory).

    For 2D embeddings (e.g. UMAP), velocity vectors are converted to
    absolute angles via atan2 and binned into [0, 2pi). For higher-
    dimensional embeddings (e.g. PCA), angles between each neighbor's
    velocity and the cell's own velocity are computed via cosine
    similarity and binned into [0, pi].

    Cells whose velocity magnitude is below 1e-10 (effectively zero)
    are assigned NaN, as their angular entropy is undefined.

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
        # 2D: absolute angles via atan2, binned into [0, 2pi)
        angles = np.arctan2(V[:, 1], V[:, 0]) % (2 * np.pi)
        bin_edges = np.linspace(0, 2 * np.pi, n_bins + 1)

        for i in tqdm(range(n_cells), desc="Scoring angular velocity entropy", unit="cell"):
            # Skip cells with zero velocity
            if zero_vel_mask[i]:
                continue
            # Only use neighbors that themselves have nonzero velocity
            nbr_idx = indices[i][~zero_vel_mask[indices[i]]]
            if len(nbr_idx) == 0:
                continue
            neighbor_angles = angles[nbr_idx]
            counts, _ = np.histogram(neighbor_angles, bins=bin_edges)
            p = counts[counts > 0] / counts.sum()
            entropy_scores[i] = -np.sum(p * np.log2(p))
    else:
        # Higher-dim: angle between each neighbor's velocity and
        # the cell's own velocity, binned into [0, pi]
        V_unit = np.zeros_like(V)
        valid = ~zero_vel_mask
        V_unit[valid] = V[valid] / vel_norms[valid, np.newaxis]
        bin_edges = np.linspace(0, np.pi, n_bins + 1)

        for i in tqdm(range(n_cells), desc="Scoring angular velocity entropy", unit="cell"):
            if zero_vel_mask[i]:
                continue
            nbr_idx = indices[i][~zero_vel_mask[indices[i]]]
            if len(nbr_idx) == 0:
                continue
            cos_sim = V_unit[nbr_idx] @ V_unit[i]
            cos_sim = np.clip(cos_sim, -1.0, 1.0)
            theta = np.arccos(cos_sim)
            counts, _ = np.histogram(theta, bins=bin_edges)
            p = counts[counts > 0] / counts.sum()
            entropy_scores[i] = -np.sum(p * np.log2(p))

    # Normalize to [0, 1] if requested
    if normalize:
        max_entropy = np.log2(n_bins)
        if max_entropy > 0:
            entropy_scores = entropy_scores / max_entropy

    adata.obs[key_added] = entropy_scores

    if not inplace:
        return adata
