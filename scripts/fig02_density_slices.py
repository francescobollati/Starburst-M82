"""
Fig. 2 -- Time-averaged density slices of the NoSources (left) and fid (right)
runs. In each panel the left half shows the hot phase and the right half the
cold phase (cold mass divided by cell volume).

Requires the snapshots of the ``fid`` and ``NoSources`` runs.

Usage
-----
    python scripts/fig02_density_slices.py [--sim-root /path/to/Simulations]
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

from common import (SERIF, UNIT_DENS, SliceGrid, averaging_snapshots,
                    run_output_dir, savefig, set_style)

N_PIXELS = 512
PLOT_RANGE = 3.5          # kpc cut from each side of the 20 kpc box
NORM_HOT = colors.LogNorm(vmin=1e-29, vmax=1e-24)
NORM_COLD = colors.LogNorm(vmin=1e-28, vmax=1e-23)


def averaged_density_maps(run, sim_root):
    """Mean hot and cold density slices over the averaging window [g cm^-3]."""
    grid, hot, cold, count = None, 0.0, 0.0, 0
    for snap in averaging_snapshots(run_output_dir(run, sim_root)):
        if grid is None:
            grid = SliceGrid(snap.boxsize, N_PIXELS)
        idx = grid.nearest_cells(snap)
        hot = hot + grid.sample(snap.density, idx) * UNIT_DENS
        cold = cold + grid.sample(snap.cold_density, idx) * UNIT_DENS
        count += 1
    return grid, hot / count, cold / count


def draw_split_panel(ax, grid, hot, cold, title):
    """Hot phase on the left half of the box, cold phase on the right half."""
    mid = np.searchsorted(grid.edges, grid.boxsize / 2)
    pc_hot = ax.pcolormesh(grid.edges[:mid + 1], grid.edges, hot[:, :mid],
                           rasterized=True, cmap="Reds", norm=NORM_HOT)
    pc_cold = ax.pcolormesh(grid.edges[mid:], grid.edges, cold[:, mid:],
                            rasterized=True, cmap="Blues", norm=NORM_COLD)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(PLOT_RANGE, grid.boxsize - PLOT_RANGE)
    ax.set_ylim(PLOT_RANGE, grid.boxsize - PLOT_RANGE)
    ax.set_xlabel("[kpc]")
    ax.set_yticklabels([])
    ax.set_title(title, fontsize=15, family=SERIF)
    return pc_hot, pc_cold


def add_colorbar(fig, mappable, x0, pos, label=None, ticks_left=False):
    cax = fig.add_axes([x0, pos.y0, 0.012, pos.height])
    cbar = fig.colorbar(mappable, cax=cax, orientation="vertical")
    if ticks_left:
        cbar.ax.yaxis.set_ticks_position("left")
        cbar.ax.yaxis.set_label_position("left")
    if label:
        cbar.set_label(label, fontsize=11, family=SERIF)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sim-root", default=None)
    args = parser.parse_args()

    set_style(ticks_in=True)
    panels = [("NoSources", "No Source Terms"), ("fid", "Source Terms")]

    fig, axes = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={"wspace": 0.2}, dpi=300)
    label = r"$\rho$  [g/cm$^3$]"
    for ax, (run, title) in zip(axes, panels):
        grid, hot, cold = averaged_density_maps(run, args.sim_root)
        pc_hot, pc_cold = draw_split_panel(ax, grid, hot, cold, title)

        # Colorbars hugging each panel: hot on the left, cold on the right.
        pos = ax.get_position()
        outer_left, outer_right = ax is axes[0], ax is axes[-1]
        add_colorbar(fig, pc_hot, pos.x0 - 0.01, pos, label if outer_left else None,
                     ticks_left=True)
        add_colorbar(fig, pc_cold, pos.x1, pos, label if outer_right else None)

        centre = 0.5 * (pos.x0 + pos.x1)
        for dx, text in ((-0.1, "Hot phase"), (0.1, "Cold phase")):
            fig.text(centre + dx, 0.8, text, ha="center", va="bottom",
                     fontsize=13, family=SERIF)

    savefig(fig, "fig02_density_slices.png")


if __name__ == "__main__":
    main()
