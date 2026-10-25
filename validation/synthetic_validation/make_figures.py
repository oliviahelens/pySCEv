"""
Regenerate the synthetic-validation figures from the same hand-built
velocity fields used in tests/test_velocity_synthetic.py.

    synthetic_scenarios.png   aligned / bifurcation / random fields + scores
    score_distributions.png   per-cell score histograms for the three fields
    mixed_dataset.png         aligned + random + zero-velocity clusters
    rotation_invariance.png   coherent flow at different headings
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from anndata import AnnData

import pysce

HERE = Path(__file__).resolve().parent
KEY = "angular_velocity_entropy"


def scored(pos, vel, **kw):
    adata = AnnData(np.zeros((len(pos), 1)))
    adata.obsm["X_umap"] = pos
    adata.obsm["velocity_umap"] = vel
    pysce.score_angular_velocity_entropy(adata, **kw)
    return adata.obs[KEY].to_numpy()


def unit(theta):
    return np.column_stack([np.cos(theta), np.sin(theta)])


def main():
    rng = np.random.default_rng(0)
    n = 400
    pos = rng.normal(0, 0.1, (n, 2))
    fields = {
        "aligned": unit(np.zeros(n)),
        "bifurcation (±45°)": unit(np.where(np.arange(n) < n // 2, np.pi / 4, -np.pi / 4)),
        "random": unit(rng.uniform(0, 2 * np.pi, n)),
    }
    scores = {k: scored(pos, v, n_neighbors=30) for k, v in fields.items()}

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (name, vel) in zip(axes, fields.items()):
        sc = ax.quiver(pos[:, 0], pos[:, 1], vel[:, 0], vel[:, 1], scores[name],
                       cmap="viridis", clim=(0, 1), scale=25, width=0.004)
        ax.set_title(f"{name}: median {np.nanmedian(scores[name]):.2f}")
        ax.set_aspect("equal")
    fig.colorbar(sc, ax=axes, label="angular velocity entropy")
    fig.savefig(HERE / "synthetic_scenarios.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    for name, s in scores.items():
        ax.hist(s, bins=np.linspace(0, 1, 41), alpha=0.6, label=name)
    ax.set_xlabel("angular velocity entropy")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(HERE / "score_distributions.png", dpi=150)
    plt.close(fig)

    m = 200
    pos_m = np.vstack([rng.normal([-5, 0], 0.3, (m, 2)), rng.normal([5, 0], 0.3, (m, 2)),
                       rng.normal([0, 5], 0.3, (m, 2))])
    vel_m = np.vstack([unit(np.zeros(m)), unit(rng.uniform(0, 2 * np.pi, m)), np.zeros((m, 2))])
    s = scored(pos_m, vel_m)
    fig, ax = plt.subplots(figsize=(6, 5))
    nan = np.isnan(s)
    ax.scatter(*pos_m[nan].T, c="lightgrey", s=8, label="zero velocity (NaN)")
    sc = ax.scatter(*pos_m[~nan].T, c=s[~nan], cmap="viridis", vmin=0, vmax=1, s=8)
    fig.colorbar(sc, label="angular velocity entropy")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("mixed: aligned (left), random (right), zero (top)")
    fig.tight_layout()
    fig.savefig(HERE / "mixed_dataset.png", dpi=150)
    plt.close(fig)

    headings = np.arange(0, 91, 2.5)
    med = [np.median(scored(pos, unit(np.deg2rad(h) + rng.normal(0, np.deg2rad(5), n))))
           for h in headings]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(headings, med, marker="o")
    ax.set_ylim(0, 0.5)
    ax.set_xlabel("flow heading (degrees)")
    ax.set_ylabel("median entropy")
    ax.set_title("Coherent flow (±5° noise): score does not depend on heading")
    fig.tight_layout()
    fig.savefig(HERE / "rotation_invariance.png", dpi=150)
    plt.close(fig)
    print("bifurcation median", np.nanmedian(scores["bifurcation (±45°)"]),
          "random median", np.nanmedian(scores["random"]), "heading range", min(med), max(med))


if __name__ == "__main__":
    main()
