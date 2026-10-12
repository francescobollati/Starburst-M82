"""
Fig. A.3 -- Ratio between the source terms of the power-law (ppl) and delta (d)
cloud models: mass loss, mass growth and drag force, as a function of the
critical radius r_crit (power law) and of the cloud size r_cl (delta).

Solid black lines mark where the two models agree (ratio = 1). In the middle
panel the dashed line repeats the mass-loss equality line of the first panel.
The crosses mark the rcl0.01 and rcl0.6 runs at the median r_crit of fid.

Usage
-----
    python scripts/figA3_powerlaw_vs_delta.py
"""

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable

from cloud_model import (mass_to_radius, normalisation, radius_to_mass, rinv_drag,
                         rinv_growth, rinv_loss)
from common import KPC, MSUN, SERIF, savefig, set_style

M_MIN, M_MAX, GAMMA = 1e-2 * MSUN, 1e9 * MSUN, 0.0
RHO_CL = 1e-25          # [g cm^-3]
N_GRID = 100
LABEL_SIZE = 14


# --- ratios S^ppl / S^d of the cell source terms (both per unit cold mass) ----
def ratio_loss(r_crit, r_cl):
    m_crit = radius_to_mass(r_crit, RHO_CL)
    return normalisation(m_crit, M_MIN, M_MAX, GAMMA) * rinv_loss(m_crit, M_MIN, M_MAX, GAMMA, RHO_CL) * r_cl


def ratio_growth(r_crit, r_cl):
    m_crit = radius_to_mass(r_crit, RHO_CL)
    xi = r_cl / r_crit
    alpha = np.where(r_cl < r_crit, 1 / 2, 1 / 4)
    return (normalisation(m_crit, M_MIN, M_MAX, GAMMA)
            * rinv_growth(m_crit, M_MIN, M_MAX, GAMMA, RHO_CL) * r_cl / xi**alpha)


def ratio_drag(r_crit, r_cl):
    m_crit = radius_to_mass(r_crit, RHO_CL)
    m_cl = radius_to_mass(r_cl, RHO_CL)
    return (normalisation(m_crit, M_MIN, M_MAX, GAMMA)
            * rinv_drag(m_crit, M_MIN, M_MAX, GAMMA, RHO_CL) * m_cl / r_cl**2)


def equality_line(ax, x, y, ratio, **kwargs):
    """Draw the ratio = 1 contour and return its vertices."""
    contour = ax.contour(x, y, ratio, levels=[1], colors="none")
    vertices = np.concatenate(contour.allsegs[0])
    ax.plot(vertices[:, 0], vertices[:, 1], color="k", lw=1, **kwargs)
    return vertices


def annotate(ax, text_ratio, symbol):
    ax.text(0.08, 0.92, text_ratio, transform=ax.transAxes, ha="left", va="top",
            color="white", fontsize=12, fontname=SERIF)
    ax.add_patch(FancyArrowPatch((3e-3, 1e-2), (3e-3, 1e-1), arrowstyle="->",
                                 color="white", lw=1.5, mutation_scale=15))
    ax.text(4e-3, 1e-2 * 10**0.5, rf"decreasing ${symbol}^{{ppl}}$", color="white",
            fontsize=LABEL_SIZE, va="center", ha="left", rotation=90, fontname=SERIF)
    ax.add_patch(FancyArrowPatch((2e-1, 3e-3), (2, 3e-3), arrowstyle="->",
                                 color="white", lw=1.5, mutation_scale=15))
    ax.text(2e-1, 5e-3, rf"decreasing ${symbol}^{{d}}$", color="white",
            fontsize=LABEL_SIZE, va="center", ha="left", fontname=SERIF)


def main():
    set_style()
    r_min = mass_to_radius(M_MIN, RHO_CL) / KPC
    radii = np.logspace(np.log10(r_min), np.log10(5), N_GRID) * KPC
    r_cl, r_crit = np.meshgrid(radii, radii)
    x, y = r_cl / KPC, r_crit / KPC

    levels = np.logspace(-4, 4, 1000)
    norm = mcolors.LogNorm(vmin=1e-3, vmax=1e3)
    fig, ax = plt.subplots(1, 3, figsize=(18, 6), gridspec_kw={"wspace": 0.2}, dpi=300)

    panels = [
        (ratio_loss, r"$\dot{M}^{ppl}_{loss}/ \dot{M}^{d}_{loss}$", r"\dot{M}"),
        (ratio_growth, r"$\dot{M}^{ppl}_{growth}/ \dot{M}^{d}_{growth}$", r"\dot{M}"),
        (ratio_drag, r"$\dot{P}^{ppl}_{drag}/ \dot{P}^{d}_{drag}$", r"\dot{P}"),
    ]
    for i, (ratio_fn, text, symbol) in enumerate(panels):
        ratio = ratio_fn(r_crit, r_cl)
        ax[i].contourf(x, y, ratio, levels=levels, cmap="seismic", norm=norm)
        if i == 1:
            ax[i].plot(loss_line[:, 0], loss_line[:, 1], color="k", lw=1, linestyle="--")
        line = equality_line(ax[i], x, y, ratio)
        if i == 0:
            loss_line = line
        ax[i].set_xscale("log")
        ax[i].set_yscale("log")
        ax[i].set_xlabel("$r_{cl}$  [kpc]")
        annotate(ax[i], text, symbol)

        top = ax[i].twiny()
        top.set_xscale("log")
        top.set_xlim(radius_to_mass(radii[0], RHO_CL) / MSUN, radius_to_mass(radii[-1], RHO_CL) / MSUN)
        top.set_xlabel(r"$m_{cl}$   [$M_\odot$]")

    ax[0].set_ylabel("$r_{crit}$  [kpc]")
    # rcl0.01 and rcl0.6 runs, at the typical r_crit of the fiducial run.
    ax[0].scatter([0.01, 0.6], [0.015, 0.015], marker="x", color="k", s=50)

    # Shared colorbar on the right of the last panel. Note: as in the
    # published figure, the bar spans 1e-4-1e4 while the filled contours
    # saturate at 1e-3 and 1e3.
    sm = cm.ScalarMappable(cmap="seismic", norm=mcolors.LogNorm(vmin=levels.min(), vmax=levels.max()))
    cax = make_axes_locatable(ax[2]).append_axes("right", size="3.5%", pad=0.02)
    fig.colorbar(sm, cax=cax).set_label(r"$\dot{S}^{ppl}/\dot{S}^{d}$")

    savefig(fig, "figA3_powerlaw_vs_delta.png")


if __name__ == "__main__":
    main()
