#!/usr/bin/env python3
"""
Analyze cluster_data.root produced by the Heed / Garfield++ simulation.

Reads the eClusterTree and produces:
  - 1D histograms of energy for clusters, electrons, ions, photons
  - 2D XY-plane distributions for clusters, electrons, ions, photons
  - 3D scatter plot of cluster positions

Outputs are saved as PNG files in ./plots/.

Notes about the source ROOT file
--------------------------------
The C++ producer has several quirks the reader must tolerate:
  * Several TBranches were declared with the wrong leaflist
    (`"cluster.x/D"` was passed for y and z), so the on-disk
    branch names may collapse or be renamed by ROOT.
  * `cluster_ion` and `cluster_photon` were both registered under
    the name `cluster_electron`, so only the first one to be
    registered survives under that name; the others may be missing
    or merged.
  * `treecluster->Fill()` is called inside every inner loop, so each
    TTree entry is "one of" {cluster summary, an ion, a photon, an
    electron}, with stale values in the other fields. This script
    treats every entry as an independent point and lets the reader
    decide which fields are meaningful per entry. In practice the
    cluster-level (x, y, z, energy) values stay valid across the
    inner loops because they were assigned at the top of the outer
    loop and not overwritten, so the per-entry positions still trace
    out the cluster locations.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

try:
    import uproot
except ImportError:  # pragma: no cover
    sys.stderr.write(
        "uproot is required. Install with: pip install uproot awkward\n"
    )
    raise


# --------------------------------------------------------------------------- #
# Configuration                                                               #
# --------------------------------------------------------------------------- #

ROOT_FILE = Path("cluster_data.root")
TREE_NAME = "eClusterTree"
OUT_DIR = Path("./png")

# Histogram binning
ENERGY_BINS = 100
XY_BINS = 100


# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #

def find_branch(tree, *candidates: str) -> str | None:
    """
    Return the first branch name from *candidates* that exists in *tree*.

    Also tries case-insensitive and dot/underscore-tolerant matches, since
    the producer uses names like 'cluster.x' and 'cluster_electron.Electron_x'
    which ROOT or uproot may surface in slightly different forms.
    """
    available = list(tree.keys())
    available_norm = {k.lower().replace(".", "_"): k for k in available}
    for c in candidates:
        if c in available:
            return c
        key = c.lower().replace(".", "_")
        if key in available_norm:
            return available_norm[key]
    return None


def read_array(tree, *candidates: str) -> np.ndarray:
    """Read the first existing branch and return its values as a flat numpy array."""
    name = find_branch(tree, *candidates)
    if name is None:
        return np.array([])
    arr = tree[name].array(library="np")
    # Some branches may come back as awkward / object arrays of arrays;
    # flatten defensively.
    try:
        return np.asarray(arr, dtype=float).ravel()
    except (TypeError, ValueError):
        flat = []
        for v in arr:
            if hasattr(v, "__iter__"):
                flat.extend(v)
            else:
                flat.append(v)
        return np.asarray(flat, dtype=float)


def safe_filter(values: np.ndarray, *companions: np.ndarray):
    """
    Drop entries where the value is exactly 0 *and* every companion is 0.

    The producer initialises all fields to 0 and only overwrites the ones
    relevant to the current inner loop, so genuine non-zero rows are the
    useful ones. We keep a row if any of (value, *companions) is non-zero,
    which keeps real 0-energy / 0-position points if they ever occur in
    combination with non-zero coordinates.
    """
    if len(values) == 0:
        return values, *companions
    stack = np.vstack([values, *companions])
    mask = np.any(stack != 0, axis=0)
    return (values[mask], *(c[mask] for c in companions))


def save_hist1d(values, title, xlabel, filename, bins=ENERGY_BINS, color="steelblue"):
    """Plot and save a 1D histogram."""
    fig, ax = plt.subplots(figsize=(8, 5))
    if len(values) == 0:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes, fontsize=14)
    else:
        ax.hist(values, bins=bins, color=color, edgecolor="black", alpha=0.85)
        ax.set_yscale("log")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = OUT_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}  (n={len(values)})")


def save_hist2d(x, y, title, filename, bins=XY_BINS, cmap="viridis"):
    """Plot and save a 2D XY-plane histogram."""
    fig, ax = plt.subplots(figsize=(7, 6))
    if len(x) == 0 or len(y) == 0:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center",
                transform=ax.transAxes, fontsize=14)
    else:
        h, xedges, yedges, im = ax.hist2d(
            x, y, bins=bins, cmap=cmap, cmin=1
        )
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Count")
    ax.set_title(title)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = OUT_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}  (n={len(x)})")


def save_scatter3d(x, y, z, title, filename, color_by=None):
    """Plot and save a 3D scatter of cluster positions."""
    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    if len(x) == 0:
        ax.text2D(0.5, 0.5, "No data available", ha="center", va="center",
                  transform=ax.transAxes, fontsize=14)
    else:
        if color_by is not None and len(color_by) == len(x):
            sc = ax.scatter(x, y, z, c=color_by, cmap="plasma",
                            s=12, alpha=0.7, depthshade=True)
            cbar = fig.colorbar(sc, ax=ax, pad=0.1, shrink=0.7)
            cbar.set_label("Energy (eV)")
        else:
            ax.scatter(x, y, z, s=12, alpha=0.7, c="crimson", depthshade=True)
    ax.set_title(title)
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel("z (cm)")
    fig.tight_layout()
    out = OUT_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}  (n={len(x)})")


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #

def main():
    if not ROOT_FILE.exists():
        sys.exit(f"ROOT file not found: {ROOT_FILE.resolve()}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Opening {ROOT_FILE} ...")
    with uproot.open(ROOT_FILE) as f:
        if TREE_NAME not in f:
            sys.exit(
                f"Tree '{TREE_NAME}' not found. Keys present: {list(f.keys())}"
            )
        tree = f[TREE_NAME]

        print(f"Tree '{TREE_NAME}' has {tree.num_entries} entries.")
        print("Available branches:")
        for k in tree.keys():
            print(f"  - {k}")
        print()

        # ----- Cluster-level scalars ---------------------------------------
        # These were declared with broken leaflists; try several names.
        cl_x = read_array(tree, "cluster.x", "cluster_x")
        cl_y = read_array(tree, "cluster.y", "cluster_y", "cluster.x_1")
        cl_z = read_array(tree, "cluster.z", "cluster_z", "cluster.x_2")
        cl_t = read_array(tree, "cluster.t", "cluster_t")
        cl_e = read_array(tree, "cluster.energy", "cluster_energy")

        # If y/z fell back to the same buffer as x (because of the leaflist
        # bug), they'll all be identical; warn the user but continue.
        if len(cl_y) and len(cl_x) and np.array_equal(cl_x, cl_y):
            print("WARNING: cluster.x and cluster.y appear identical — "
                  "likely caused by the duplicated leaflist 'cluster.x/D' "
                  "in the C++ producer.")

        # ----- Electron / Ion / Photon sub-records -------------------------
        # The C++ side registered all three struct branches under the name
        # 'cluster_electron', so on disk we may see either:
        #   * one branch 'cluster_electron' (last write wins), or
        #   * three branches with suffixes added by ROOT.
        # We try both naming conventions.
        e_x = read_array(tree, "cluster_electron/Electron_x",
                         "cluster_electron.Electron_x")
        e_y = read_array(tree, "cluster_electron/Electron_y",
                         "cluster_electron.Electron_y")
        e_z = read_array(tree, "cluster_electron/Electron_z",
                         "cluster_electron.Electron_z")
        e_e = read_array(tree, "cluster_electron/Electron_energy",
                         "cluster_electron.Electron_energy")

        i_x = read_array(tree, "cluster_ion/Ion_x", "cluster_ion.Ion_x")
        i_y = read_array(tree, "cluster_ion/Ion_y", "cluster_ion.Ion_y")
        i_z = read_array(tree, "cluster_ion/Ion_z", "cluster_ion.Ion_z")
        i_e = read_array(tree, "cluster_ion/Ion_energy", "cluster_ion.Ion_energy")

        p_x = read_array(tree, "cluster_photon/Photon_x",
                         "cluster_photon.Photon_x")
        p_y = read_array(tree, "cluster_photon/Photon_y",
                         "cluster_photon.Photon_y")
        p_z = read_array(tree, "cluster_photon/Photon_z",
                         "cluster_photon.Photon_z")
        p_e = read_array(tree, "cluster_photon/Photon_energy",
                         "cluster_photon.Photon_energy")

    # Drop the "all-zero placeholder" rows that come from the producer's
    # pattern of filling the tree inside every inner loop.
    cl_e_f, cl_x_f, cl_y_f, cl_z_f = safe_filter(cl_e, cl_x, cl_y, cl_z)
    # Time is filtered on its own — use energy as the "is this a real
    # cluster row" companion so that t == 0 entries are kept only when
    # they coincide with a non-zero energy (which would be a genuine
    # zero-time cluster, e.g. the primary at t=0).
    if len(cl_t) and len(cl_e):
        cl_t_f, _ = safe_filter(cl_t, cl_e)
    else:
        cl_t_f = cl_t
    e_e_f, e_x_f, e_y_f, e_z_f = safe_filter(e_e, e_x, e_y, e_z)
    i_e_f, i_x_f, i_y_f, i_z_f = safe_filter(i_e, i_x, i_y, i_z)
    p_e_f, p_x_f, p_y_f, p_z_f = safe_filter(p_e, p_x, p_y, p_z)

    # ------------------------------------------------------------------ #
    # 1D energy histograms                                               #
    # ------------------------------------------------------------------ #
    print("\nGenerating 1D energy histograms ...")
    save_hist1d(cl_e_f, "Cluster Energy Distribution",
                "Energy (eV)", "energy_cluster.png", color="steelblue")
    save_hist1d(e_e_f, "Secondary Electron Energy Distribution",
                "Energy (eV)", "energy_electron.png", color="seagreen")
    save_hist1d(i_e_f, "Ion Energy Distribution",
                "Energy (eV)", "energy_ion.png", color="indianred")
    save_hist1d(p_e_f, "Photon Energy Distribution",
                "Energy (eV)", "energy_photon.png", color="darkorange")

    # ------------------------------------------------------------------ #
    # Cluster time distribution                                          #
    # ------------------------------------------------------------------ #
    print("\nGenerating cluster time distribution ...")
    save_hist1d(cl_t_f, "Cluster Time Distribution",
                "Time (ns)", "time_cluster.png",
                bins=ENERGY_BINS, color="slateblue")

    # ------------------------------------------------------------------ #
    # 2D XY distributions                                                #
    # ------------------------------------------------------------------ #
    print("\nGenerating 2D XY distributions ...")
    save_hist2d(cl_x_f, cl_y_f, "Cluster XY Distribution",
                "xy_cluster.png", cmap="viridis")
    save_hist2d(e_x_f, e_y_f, "Electron XY Distribution",
                "xy_electron.png", cmap="Greens")
    save_hist2d(i_x_f, i_y_f, "Ion XY Distribution",
                "xy_ion.png", cmap="Reds")
    save_hist2d(p_x_f, p_y_f, "Photon XY Distribution",
                "xy_photon.png", cmap="Oranges")

    # ------------------------------------------------------------------ #
    # 3D scatter of cluster positions                                    #
    # ------------------------------------------------------------------ #
    print("\nGenerating 3D cluster scatter ...")
    save_scatter3d(cl_x_f, cl_y_f, cl_z_f,
                   "Cluster Positions in 3D",
                   "cluster_3d.png",
                   color_by=cl_e_f if len(cl_e_f) == len(cl_x_f) else None)

    print(f"\nDone. All plots saved to {OUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
