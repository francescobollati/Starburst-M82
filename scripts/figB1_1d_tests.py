"""
Fig. B.1 -- One-dimensional stationary multiphase outflow test (Appendix B):
densities, velocities and hot-phase pressure from the AREPO test runs
(markers) compared with the semi-analytic solutions (lines), for the delta
and power-law cloud models.

The semi-analytic solutions were obtained with the code of Fielding & Bryan
(2022), modified to include our cloud models. Data in ``data/1d_tests/``; each
array has rows [R (kpc), rho_hot, rho_cold (g cm^-3), v_hot, v_cold (km/s),
P_hot (bar)].

Usage
-----
    python scripts/figB1_1d_tests.py
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D

from common import DATA_DIR, savefig, set_style

STEP = 10          # plot every 10th simulation cell
MARKER_SIZE = 30
HOT, COLD = ("r", "lightsalmon"), ("b", "turquoise")   # (line colour, marker colour)


def load(name):
    return np.load(os.path.join(DATA_DIR, "1d_tests", name + ".npy"))


def plot_quantity(ax, row, line_color, marker_color, delta_th, delta_sim, ppl_th, ppl_sim):
    ax.plot(delta_th[0], delta_th[row], color=line_color)
    ax.scatter(delta_sim[0, ::STEP], delta_sim[row, ::STEP], color=marker_color,
               marker="x", s=MARKER_SIZE)
    ax.plot(ppl_th[0], ppl_th[row], color=line_color, linestyle=":")
    ax.scatter(ppl_sim[0, ::STEP], ppl_sim[row, ::STEP], facecolors="none",
               edgecolors=marker_color, marker="o", s=MARKER_SIZE)


def handle(color, linestyle="-", marker=None, filled=True):
    if marker is None:
        return Line2D([0], [0], color=color, linestyle=linestyle)
    return Line2D([0], [0], linestyle="None", marker=marker, markeredgecolor=color,
                  markerfacecolor=color if filled else "none", color=color)


def phase_handles(line_color, marker_color):
    """Legend entry showing all four line/marker styles of one phase."""
    return (handle(line_color), handle(marker_color, marker="x"),
            handle(line_color, ":"), handle(marker_color, marker="o", filled=False))


def combined_legend(ax, groups, labels, loc):
    ax.legend(groups, labels, handler_map={tuple: HandlerTuple(ndivide=None)},
              loc=loc, fontsize=10)


def main():
    set_style()
    data = [load(n) for n in ("delta_teo", "delta_sim", "ppl_teo", "ppl_sim")]

    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(6.9, 6.0), dpi=300)
    fig.subplots_adjust(left=0.13, bottom=0.09, right=0.98, top=0.98)

    # Density
    plot_quantity(ax[0], 1, *HOT, *data)
    plot_quantity(ax[0], 2, *COLD, *data)
    combined_legend(ax[0], [phase_handles(*HOT), phase_handles(*COLD)],
                    ["Hot", "Cold"], "upper right")
    ax[0].set_ylabel(r"$\rho$ [g cm$^{-3}$]")
    ax[0].set_yscale("log")
    ax[0].set_ylim(1e-29, 1e-25)

    # Velocity
    plot_quantity(ax[1], 3, *HOT, *data)
    plot_quantity(ax[1], 4, *COLD, *data)
    combined_legend(ax[1], [(handle("k"), handle("k", ":")),
                            (handle("k", marker="x"), handle("k", marker="o", filled=False))],
                    ["semi-analytic", "simulation"], "lower right")
    ax[1].set_ylabel(r"$v$ [km s$^{-1}$]")
    ax[1].set_yscale("log")
    ax[1].set_ylim(100, 5000)

    # Hot-phase pressure
    plot_quantity(ax[2], 5, *HOT, *data)
    combined_legend(ax[2], [(handle("k"), handle("k", marker="x")),
                            (handle("k", ":"), handle("k", marker="o", filled=False))],
                    ["delta", "power law"], "upper right")
    ax[2].set_ylabel(r"$P_{hot}$   [bar]")
    ax[2].set_yscale("log")
    ax[2].set_ylim(1e-14, 1e-10)
    ax[2].set_xlabel(r"$R$   [kpc]")
    ax[2].set_xlim(0.0, 6.2)

    savefig(fig, "figB1_1d_tests.png")


if __name__ == "__main__":
    main()
