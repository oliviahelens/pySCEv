"""
Primary pancreas pass for angular velocity entropy, plus the confound
controls that the first pass lacked.

Produces every figure in this folder (except alpha_phase_portraits.png,
which checks velocity rather than the entropy metric) and a summary
file with the numbers quoted in README.md.

Controls
--------
1. Speed. Slow cells have poorly determined directions, so their angles
   are mostly noise. We (a) rescore with min_speed set to the 25th
   percentile of speed and (b) regress entropy on log speed and look at
   the per-celltype residuals. If the celltype ranking survives both,
   it is not just a speed ranking.
2. Permutation. Shuffle velocity vectors across cells (N times) and
   rescore. Measures how much celltype separation (Kruskal-Wallis H) and
   spatial smoothness (Moran's I on the UMAP k-NN graph) the scoring
   kernel produces on its own.

Usage
-----
    python run_analysis.py                      # downloads data via scVelo
    python run_analysis.py --data cache/endocrinogenesis_day15.h5ad
    python run_analysis.py --permute 100
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import pysce  # noqa: E402
from _common.score_pipeline import ScoreParams, score  # noqa: E402

KEY = "angular_velocity_entropy"
ORDER = ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine",
         "Beta", "Alpha", "Delta", "Epsilon"]


# --------------------------------------------------------------------------
# Data + scoring
# --------------------------------------------------------------------------
def load(path: Path | None):
    import scvelo as scv
    cache = HERE / "cache" / "endocrinogenesis_day15.h5ad"
    path = path or cache
    adata = scv.datasets.pancreas(file_path=str(path))
    print(f"[load] {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def expression_entropy(adata, cache: Path) -> np.ndarray:
    """pySCE expression entropy on the full count matrix (cached; slow on CPU)."""
    if cache.exists():
        return np.load(cache)
    from anndata import AnnData
    full = AnnData(adata.layers["spliced"] + adata.layers["unspliced"])
    full.var_names = adata.var_names.str.upper()
    full.var_names_make_unique()
    pysce.score(full, use_raw=False, batch_size=2)
    ent = full.obs["entropy"].to_numpy()
    cache.parent.mkdir(exist_ok=True)
    np.save(cache, ent)
    return ent


def knn(X, k):
    from sklearn.neighbors import NearestNeighbors
    return NearestNeighbors(n_neighbors=k).fit(X).kneighbors(return_distance=False)


def mean_angular_deviation(V, idx, min_speed=1e-10):
    """Mean angle between a cell's velocity and each neighbor's (first moment)."""
    norms = np.linalg.norm(V, axis=1)
    valid = norms >= min_speed
    U = np.zeros_like(V)
    U[valid] = V[valid] / norms[valid, None]
    cos = np.clip(np.einsum("nd,nkd->nk", U, U[idx]), -1, 1)
    ang = np.where(valid[idx], np.arccos(cos), np.nan)
    out = np.nanmean(ang, axis=1)
    out[~valid] = np.nan
    return out


def morans_i(x, idx):
    ok = np.isfinite(x)
    z = np.where(ok, x - np.nanmean(x), 0.0)
    w = ok[:, None] & ok[idx]
    num = np.sum(z[:, None] * z[idx] * w) / w.sum()
    return float(num / np.mean(z[ok] ** 2))


def kruskal(x, groups):
    ok = np.isfinite(x)
    samples = [x[ok & (groups == g)] for g in ORDER]
    return float(stats.kruskal(*[s for s in samples if len(s) >= 15]).statistic)


def medians(x, groups):
    return {g: float(np.nanmedian(x[groups == g])) for g in ORDER}


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
def palette(adata):
    cols = adata.uns.get("clusters_colors")
    cats = list(adata.obs["clusters"].cat.categories)
    if cols is not None and len(cols) == len(cats):
        return dict(zip(cats, cols))
    cmap = plt.get_cmap("tab10")
    return {c: cmap(i) for i, c in enumerate(cats)}


def scatter_by_type(ax, adata, x, y, pal):
    ct = adata.obs["clusters"].astype(str).to_numpy()
    for g in ORDER:
        m = (ct == g) & np.isfinite(x) & np.isfinite(y)
        ax.scatter(x[m], y[m], s=6, alpha=0.6, color=pal[g], label=g)


def box_by_type(ax, adata, y, title):
    ct = adata.obs["clusters"].astype(str).to_numpy()
    data = [y[(ct == g) & np.isfinite(y)] for g in ORDER]
    ax.boxplot(data, showfliers=False)
    ax.set_xticks(range(1, len(ORDER) + 1), ORDER, rotation=40, ha="right")
    ax.set_title(title)


def fig_umap_panels(adata, out):
    import scvelo as scv
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    scv.pl.scatter(adata, color="clusters", ax=axes[0], show=False,
                   legend_loc="right margin", title="cell type")
    scv.pl.scatter(adata, color=KEY, ax=axes[1], show=False,
                   color_map="viridis", title="angular velocity entropy")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_streamlines(adata, out):
    import scvelo as scv
    fig, ax = plt.subplots(figsize=(7, 6))
    scv.pl.velocity_embedding_stream(
        adata, basis="umap", color=KEY, color_map="viridis", ax=ax,
        show=False, title="velocity streamlines over angular velocity entropy",
    )
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_box(adata, y, out, title):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    box_by_type(ax, adata, y, title)
    ax.set_ylabel("angular velocity entropy")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_scatter(adata, x, y, xlabel, ylabel, title, out, pal):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    scatter_by_type(ax, adata, x, y, pal)
    ok = np.isfinite(x) & np.isfinite(y)
    r = stats.pearsonr(x[ok], y[ok])[0]
    rho = stats.spearmanr(x[ok], y[ok])[0]
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\nPearson r = {r:.2f}, Spearman rho = {rho:.2f}")
    ax.legend(fontsize=7, frameon=False, markerscale=2)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_expression_umaps(adata, out):
    import scvelo as scv
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    scv.pl.scatter(adata, color="expression_entropy", ax=axes[0], show=False,
                   color_map="magma", title="expression entropy (pySCE)")
    scv.pl.scatter(adata, color=KEY, ax=axes[1], show=False,
                   color_map="viridis", title="angular velocity entropy")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_speed_control(adata, speed, resid, filtered, out, pal):
    ct = adata.obs["clusters"].astype(str).to_numpy()
    ent = adata.obs[KEY].to_numpy()
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))

    ok = np.isfinite(ent) & (speed > 0)
    scatter_by_type(axes[0], adata, np.log10(np.where(speed > 0, speed, np.nan)), ent, pal)
    rho = stats.spearmanr(np.log10(speed[ok]), ent[ok])[0]
    axes[0].set_xlabel("log10 velocity speed (UMAP)")
    axes[0].set_ylabel("angular velocity entropy")
    axes[0].set_title(f"entropy vs speed (Spearman rho = {rho:.2f})")
    axes[0].legend(fontsize=7, frameon=False, markerscale=2)

    x = np.arange(len(ORDER))
    for vals, lab, mk in [(ent, "all cells", "o"),
                          (filtered, "min_speed = 25th pct", "s")]:
        axes[1].plot(x, [np.nanmedian(vals[ct == g]) for g in ORDER],
                     marker=mk, label=lab)
    axes[1].set_xticks(x, ORDER, rotation=40, ha="right")
    axes[1].set_ylabel("median entropy")
    axes[1].set_title("celltype medians with slow cells removed")
    axes[1].legend(frameon=False)

    box_by_type(axes[2], adata, resid, "entropy residual after regressing on log speed")
    axes[2].axhline(0, color="grey", lw=0.8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig_permutation(null_kw, obs_kw, null_I, obs_I, out):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, null, obs, lab in [(axes[0], null_kw, obs_kw, "Kruskal-Wallis H (celltype)"),
                               (axes[1], null_I, obs_I, "Moran's I (UMAP k-NN)")]:
        ax.hist(null, bins=25, color="grey", alpha=0.7, label="shuffled velocity")
        ax.axvline(obs, color="crimson", lw=2, label=f"observed = {obs:.2f}")
        z = (obs - null.mean()) / null.std()
        ax.set_title(f"{lab}\nnull {null.mean():.2f} ± {null.std():.2f}, z = {z:.1f}")
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=None)
    p.add_argument("--permute", type=int, default=100)
    p.add_argument("--skip-expression-entropy", action="store_true")
    args = p.parse_args()

    adata = load(args.data)
    params = ScoreParams()
    ee = None
    if not args.skip_expression_entropy:
        print("[score] expression entropy (pySCE)")
        ee = expression_entropy(adata, HERE / "cache" / "expression_entropy.npy")

    print("[score] scVelo + angular velocity entropy")
    score(adata, params)
    if ee is not None:
        adata.obs["expression_entropy"] = ee

    pal = palette(adata)
    ct = adata.obs["clusters"].astype(str).to_numpy()
    ent = adata.obs[KEY].to_numpy()
    V = np.asarray(adata.obsm["velocity_umap"])
    speed = np.linalg.norm(V, axis=1)
    idx = knn(adata.obsm["X_umap"], params.n_neighbors)
    mad = mean_angular_deviation(V, idx)
    conf = adata.obs["velocity_confidence"].to_numpy()

    summary = {"n_cells": int(adata.n_obs),
               "n_scored": int(np.isfinite(ent).sum()),
               "celltype_median": medians(ent, ct)}

    def corr(a, b):
        ok = np.isfinite(a) & np.isfinite(b)
        return {"pearson": float(stats.pearsonr(a[ok], b[ok])[0]),
                "spearman": float(stats.spearmanr(a[ok], b[ok])[0])}

    summary["vs_velocity_confidence"] = corr(ent, conf)
    summary["vs_mean_angular_deviation"] = corr(ent, mad)
    summary["vs_log_speed"] = corr(ent, np.log10(np.where(speed > 0, speed, np.nan)))
    if ee is not None:
        summary["vs_expression_entropy"] = corr(ent, ee)

    # --- speed control -----------------------------------------------------
    q25 = float(np.quantile(speed[speed > 0], 0.25))
    pysce.score_angular_velocity_entropy(
        adata, basis="umap", n_neighbors=params.n_neighbors,
        n_bins=params.n_bins, min_speed=q25, key_added="ave_fast",
    )
    filtered = adata.obs["ave_fast"].to_numpy()
    ok = np.isfinite(ent) & (speed > 0)
    slope, icpt = np.polyfit(np.log10(speed[ok]), ent[ok], 1)
    resid = np.full_like(ent, np.nan)
    resid[ok] = ent[ok] - (slope * np.log10(speed[ok]) + icpt)

    summary["speed_control"] = {
        "min_speed_q25": q25,
        "celltype_median_min_speed_q25": medians(filtered, ct),
        "celltype_median_speed_residual": medians(resid, ct),
        "kruskal_H_raw": kruskal(ent, ct),
        "kruskal_H_min_speed_q25": kruskal(filtered, ct),
        "kruskal_H_speed_residual": kruskal(resid, ct),
        "rank_corr_celltype_medians_raw_vs_residual": float(stats.spearmanr(
            list(medians(ent, ct).values()), list(medians(resid, ct).values()))[0]),
        "rank_corr_celltype_medians_raw_vs_min_speed": float(stats.spearmanr(
            list(medians(ent, ct).values()), list(medians(filtered, ct).values()))[0]),
    }

    # --- permutation control ----------------------------------------------
    obs_kw, obs_I = kruskal(ent, ct), morans_i(ent, idx)
    null_kw, null_I = np.empty(args.permute), np.empty(args.permute)
    null_med = np.empty((args.permute, len(ORDER)))
    for i in range(args.permute):
        perm = np.random.default_rng(i).permutation(adata.n_obs)
        adata.obsm["velocity_perm"] = V[perm]
        adata.obsm["X_perm"] = adata.obsm["X_umap"]
        pysce.score_angular_velocity_entropy(
            adata, basis="perm", n_neighbors=params.n_neighbors,
            n_bins=params.n_bins, key_added="_perm")
        e = adata.obs["_perm"].to_numpy()
        null_kw[i], null_I[i] = kruskal(e, ct), morans_i(e, idx)
        null_med[i] = list(medians(e, ct).values())
    for k in ["velocity_perm", "X_perm"]:
        del adata.obsm[k]
    del adata.obs["_perm"]
    adata.uns.pop("_perm_params", None)

    summary["permutation"] = {
        "n": args.permute,
        "kruskal_H": {"observed": obs_kw, "null_mean": float(null_kw.mean()),
                      "null_sd": float(null_kw.std()),
                      "z": float((obs_kw - null_kw.mean()) / null_kw.std())},
        "morans_I": {"observed": obs_I, "null_mean": float(null_I.mean()),
                     "null_sd": float(null_I.std()),
                     "z": float((obs_I - null_I.mean()) / null_I.std())},
        "celltype_median_null_mean": dict(zip(ORDER, null_med.mean(0).round(4).tolist())),
    }

    # --- figures ------------------------------------------------------------
    fig_umap_panels(adata, HERE / "umap_panels.png")
    fig_streamlines(adata, HERE / "velocity_streamlines.png")
    fig_box(adata, ent, HERE / "entropy_by_celltype.png",
            "angular velocity entropy by cell type")
    fig_scatter(adata, conf, ent, "scVelo velocity confidence",
                "angular velocity entropy", "entropy vs velocity confidence",
                HERE / "entropy_vs_confidence.png", pal)
    fig_scatter(adata, mad, ent, "mean angular deviation (rad)",
                "angular velocity entropy", "entropy vs mean angular deviation",
                HERE / "entropy_vs_angular_deviation.png", pal)
    if ee is not None:
        fig_expression_umaps(adata, HERE / "expression_vs_velocity_entropy.png")
        fig_scatter(adata, ee, ent, "expression entropy (pySCE)",
                    "angular velocity entropy", "expression vs velocity entropy",
                    HERE / "expression_vs_velocity_scatter.png", pal)
    fig_speed_control(adata, speed, resid, filtered, HERE / "speed_control.png", pal)
    fig_permutation(null_kw, obs_kw, null_I, obs_I, HERE / "permutation_null.png")

    out = HERE / "summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    # Hand the scored object to the enrichment analysis
    cache = HERE.parent / "pancreas_endocrinogenesis_enrichment" / "cache"
    cache.mkdir(exist_ok=True)
    adata.write_h5ad(cache / "pancreas.scored.h5ad")


if __name__ == "__main__":
    main()
