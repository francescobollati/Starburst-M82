"""
Fig. A.1 -- The coefficient r_drag^-1(m_crit) of the power-law cloud model for
varying gamma, m_min and m_max (other parameters at their fiducial values),
plus a histogram of the m_crit values found in the fid simulation.

The histogram is computed from the fid snapshots with ``--sim-root`` and cached
in ``data/fid_mcrit_histogram.npz``; without snapshots the cached file is used.

Usage
-----
    python scripts/figA1_rinv_coefficients.py [--sim-root /path/to/Simulations]
"""

import argparse
import os

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from cloud_model import rinv_drag, rinv_growth, rinv_loss
from common import (DATA_DIR, KPC, MSUN, SERIF, SliceGrid, averaging_snapshots,
                    run_output_dir, savefig, set_style)

HIST_CACHE = os.path.join(DATA_DIR, "fid_mcrit_histogram.npz")
N_PIXELS = 256
CODE_MASS_IN_MSUN = 1e10

# Fiducial model parameters (Table 2) and the ranges explored in each panel.
N_CURVES = 100
M_MIN, M_MAX, GAMMA = 1e-2 * MSUN, 1e9 * MSUN, 0.0
RHO_CL = 1e-23                                       # [g cm^-3]
GAMMAS = np.linspace(-1.99, 1.99, N_CURVES)
M_MINS = np.logspace(-4, 0, N_CURVES) * MSUN
M_MAXS = np.logspace(7, 11, N_CURVES) * MSUN
M_CRIT = np.logspace(-8, 14, 1000) * MSUN

XLIM, YLIM = (1e-6, 1e14), (1e-2, 5e4)
CMAP = plt.get_cmap("copper")
CURVE_COLORS = [CMAP(i / (N_CURVES - 1)) for i in range(N_CURVES)]


def mcrit_histogram(sim_root):
    """Normalised histogram of the time-averaged m_crit slice of the fid run.

    Masses are binned in code units (1e10 Msun); the density is normalised so
    that its integral over the bins is unity.
    """
    grid, total, count = None, 0.0, 0
    for snap in averaging_snapshots(run_output_dir("fid", sim_root)):
        if grid is None:
            grid = SliceGrid(snap.boxsize, N_PIXELS)
        total = total + grid.sample(snap.cloud_mass, grid.nearest_cells(snap))
        count += 1
    values = (total / count)[grid.radius < 8.5]
    values = values[np.isfinite(values) & (values > 0)]
    counts, edges = np.histogram(values, bins=np.logspace(-16, 4, 50))
    density = counts / np.sum(counts * np.diff(edges))
    return edges, density


def load_histogram(sim_root):
    if sim_root is not None:
        edges, density = mcrit_histogram(sim_root)
        np.savez(HIST_CACHE, edges=edges, density=density)
    else:
        cached = np.load(HIST_CACHE)
        edges, density = cached["edges"], cached["density"]
    return edges, density


def add_colorbar(fig, ax, norm, label):
    sm = cm.ScalarMappable(cmap=CMAP, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax).set_label(label)


def drag_curves(ax, param_values, make_args):
    """One r_drag^-1(m_crit) curve per parameter value [kpc^-1]."""
    for color, value in zip(CURVE_COLORS, param_values):
        m_min, m_max, gamma = make_args(value)
        y = rinv_drag(M_CRIT, m_min, m_max, gamma, RHO_CL, clip=True) * KPC * RHO_CL
        ax.plot(M_CRIT / MSUN, y, color=color)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sim-root", default=None,
                        help="recompute the m_crit histogram from the fid snapshots")
    args = parser.parse_args()

    set_style()
    fig, axs = plt.subplots(4, 1, figsize=(6, 10), dpi=300)

    # (1) Histogram of m_crit in the fiducial simulation.
    edges, density = load_histogram(args.sim_root)
    axs[0].bar(edges[:-1] * CODE_MASS_IN_MSUN, density,
               width=np.diff(edges) * CODE_MASS_IN_MSUN, align="edge", alpha=0.7,
               color="LightSeaGreen")
    axs[0].set_xscale("log")
    axs[0].yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    axs[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axs[0].text(0.8e-5, density.max() * 0.8, "Normalised counts \n in fid. simulation",
                fontsize=12, color="black", fontname=SERIF)
    # Invisible colorbar, so that the histogram has the same width as the other panels.
    dummy = cm.ScalarMappable(cmap="viridis", norm=mcolors.Normalize(0, 1))
    fig.colorbar(dummy, ax=axs[0]).ax.set_visible(False)

    # (2) Varying the low-mass slope gamma.
    drag_curves(axs[1], GAMMAS, lambda g: (M_MIN, M_MAX, g))
    add_colorbar(fig, axs[1], mcolors.Normalize(GAMMAS.min(), GAMMAS.max()), r"$\gamma$")

    # (3) Varying m_min; r_loss^-1 and r_growth^-1 for the fiducial parameters.
    drag_curves(axs[2], M_MINS, lambda m: (m, M_MAX, GAMMA))
    axs[2].plot(M_CRIT / MSUN, rinv_loss(M_CRIT, M_MIN, M_MAX, GAMMA, RHO_CL, clip=True) * KPC,
                color="k", linestyle="--", linewidth=1.5, label="$r^{-1}_{loss}$")
    axs[2].plot(M_CRIT / MSUN, rinv_growth(M_CRIT, M_MIN, M_MAX, GAMMA, RHO_CL, clip=True) * KPC,
                color="k", linestyle=":", linewidth=1.5, label="$r^{-1}_{growth}$")
    axs[2].legend()
    add_colorbar(fig, axs[2], mcolors.LogNorm(M_MINS.min() / MSUN, M_MINS.max() / MSUN),
                 r"$m_{min}$  $[M_\odot]$")

    # (4) Varying m_max.
    drag_curves(axs[3], M_MAXS, lambda m: (M_MIN, m, GAMMA))
    add_colorbar(fig, axs[3], mcolors.LogNorm(M_MAXS.min() / MSUN, M_MAXS.max() / MSUN),
                 r"$m_{max}$  $[M_\odot]$")
    axs[3].set_xlabel(r"$m_{crit}$   [$M_\odot$]")

    axs[0].set_xlim(*XLIM)
    for ax in axs[1:]:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*XLIM)
        ax.set_ylim(*YLIM)
        ax.set_ylabel("$r^{-1}_{drag}$   [kpc$^{-1}$]")

    plt.tight_layout()
    savefig(fig, "figA1_rinv_coefficients.png")


if __name__ == "__main__":
    main()
