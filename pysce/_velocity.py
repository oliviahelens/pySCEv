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
    Modifies adata in-place.

    Params
    -------
    adata
        Annotated data matrix with spliced/unspliced layers.
    mode
        scVelo velocity mode: 'stochastic', 'deterministic', or 'dynamical'.
    """
    import scvelo as scv

    if 'velocity' in adata.layers:
        return

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

    Params
    -------
    adata
        Annotated data matrix. Must have velocity projection in
        adata.obsm[f'velocity_{basis}']. If missing, scVelo is
        used to compute it.
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

    velocity_key = f'velocity_{basis}'
    embedding_key = f'X_{basis}'

    # Validate embedding exists
    if embedding_key not in adata.obsm:
        raise ValueError(
            f"Embedding '{embedding_key}' not found in adata.obsm. "
            f"Run the corresponding embedding method first."
        )

    # Ensure velocity projection exists
    if velocity_key not in adata.obsm:
        try:
            import scvelo as scv
            scv.tl.velocity_embedding(adata, basis=basis)
        except ImportError:
            raise ImportError(
                f"Velocity embedding '{velocity_key}' not found and scvelo "
                "is not installed. Either precompute velocity embeddings or "
                "install scvelo: pip install scvelo"
            )

    V = np.asarray(adata.obsm[velocity_key])
    X_emb = np.asarray(adata.obsm[embedding_key])
    n_cells = adata.n_obs
    ndim = V.shape[1]

    # Build k-NN index in embedding space
    nn = NearestNeighbors(n_neighbors=min(n_neighbors + 1, n_cells), algorithm='auto')
    nn.fit(X_emb)
    indices = nn.kneighbors(return_distance=False)[:, 1:]  # exclude self

    entropy_scores = np.zeros(n_cells)

    if ndim == 2:
        # 2D: absolute angles via atan2, binned into [0, 2pi)
        angles = np.arctan2(V[:, 1], V[:, 0]) % (2 * np.pi)
        bin_edges = np.linspace(0, 2 * np.pi, n_bins + 1)

        for i in range(n_cells):
            neighbor_angles = angles[indices[i]]
            counts, _ = np.histogram(neighbor_angles, bins=bin_edges)
            total = counts.sum()
            if total > 0:
                p = counts[counts > 0] / total
                entropy_scores[i] = -np.sum(p * np.log2(p))
    else:
        # Higher-dim: angle between each neighbor's velocity and
        # the cell's own velocity, binned into [0, pi]
        norms = np.linalg.norm(V, axis=1, keepdims=True)
        norms = np.clip(norms, 1e-10, None)
        V_unit = V / norms
        bin_edges = np.linspace(0, np.pi, n_bins + 1)

        for i in range(n_cells):
            cos_sim = V_unit[indices[i]] @ V_unit[i]
            cos_sim = np.clip(cos_sim, -1.0, 1.0)
            theta = np.arccos(cos_sim)
            counts, _ = np.histogram(theta, bins=bin_edges)
            total = counts.sum()
            if total > 0:
                p = counts[counts > 0] / total
                entropy_scores[i] = -np.sum(p * np.log2(p))

    # Normalize to [0, 1] if requested
    if normalize:
        max_entropy = np.log2(n_bins)
        if max_entropy > 0:
            entropy_scores = entropy_scores / max_entropy

    if not inplace:
        adata = adata.copy()

    adata.obs[key_added] = entropy_scores

    if not inplace:
        return adata
