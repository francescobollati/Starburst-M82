"""
Figs. 3, 4, 5 -- Radial profiles of velocity, density and mass outflow rate of
the hot and cold phases for the fid, NoDrag_fmix0.03 and NoSources runs,
compared with the high-resolution simulation of Schneider et al. (2020, S20).

* velocity and density: median (and 16th-84th percentiles for fid) within
  0.3 kpc radial bins inside a 30 deg cone, from time-averaged slices;
* mass outflow rate: sum of m v / dr over all cells in the shell and cone,
  averaged over snapshots (shaded band: standard deviation in time).

The cold-phase density is the density of the individual sub-grid clouds,
rho_cl = chi * rho_hot (Eq. 11). S20 curves are digitised from their Figs. 6-8
and stored in ``data/schneider2020/``.

Requires the snapshots of the ``fid``, ``NoDrag_fmix0.03`` and ``NoSources`` runs.

Usage
-----
    python scripts/fig03_04_05_schneider_comparison.py [--sim-root /path/to/Simulations]
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, SymmetricalLogLocator

from common import (M_P, MU, SERIF, UNIT_DENS, SliceGrid, averaging_snapshots,
                    outflow_rate, radial_profile, read_digitized, run_output_dir,
                    savefig, set_style)

N_PIXELS = 256
RMAX_BINS = 8.75     # radial binning extent [kpc]
RMAX_PLOT = 5.7      # outer radius shown in the figures [kpc]
LEGEND_FONT = {"family": SERIF}


def time_averaged_run(run, sim_root):
    """Radial profiles of one run, averaged over the T = (27, 30) Myr window."""
    grid, count = None, 0
    maps = {"v_hot": 0.0, "v_cold": 0.0, "rho_hot": 0.0, "rho_cloud": 0.0}
    mdot_hot, mdot_cold = [], []

    for snap in averaging_snapshots(run_output_dir(run, sim_root)):
        if grid is None:
            grid = SliceGrid(snap.boxsize, N_PIXELS)
        idx = grid.nearest_cells(snap)
        maps["v_hot"] += grid.sample(snap.vrad, idx)
        maps["v_cold"] += grid.sample(snap.cold_vrad, idx)
        maps["rho_hot"] += grid.sample(snap.density, idx) * UNIT_DENS
        maps["rho_cloud"] += grid.sample(snap.cloud_density, idx) * UNIT_DENS

        r_mdot, mh = outflow_rate(snap, snap.mass, snap.vrad, RMAX_BINS)
        _, mc = outflow_rate(snap, snap.cold_mass, snap.cold_vrad, RMAX_BINS)
        mdot_hot.append(mh)
        mdot_cold.append(mc)
        count += 1

    prof = {}
    for key, image in maps.items():
        if key.startswith("rho"):
            image = image / (MU * M_P)            # number density [cm^-3]
            key = "n" + key[3:]
        r, median, p16, p84 = radial_profile(grid, image / count, RMAX_BINS)
        prof[key] = (median, p16, p84)
    prof["r"] = r
    prof["r_mdot"] = r_mdot
    prof["mdot_hot"] = (np.mean(mdot_hot, axis=0), np.std(mdot_hot, axis=0))
    prof["mdot_cold"] = (np.mean(mdot_cold, axis=0), np.std(mdot_cold, axis=0))
    return prof


def inner(r, *arrays):
    """Restrict arrays to the plotted radial range."""
    sel = r < RMAX_PLOT
    return [r[sel]] + [a[sel] for a in arrays]


def plot_schneider(ax, hot_name, cold_name):
    for name, color in ((hot_name, "r"), (cold_name, "b")):
        r, y = read_digitized(name)
        ax.plot(*inner(r, y), color=color, linestyle="--")


def legend_runs(include_nosources_cold=True):
    handles = [
        Line2D([], [], color="r", linestyle="-", linewidth=1.5, label="fid"),
        Line2D([], [], color="r", linestyle=":", linewidth=1.5, label="NoDrag_fmix0.03"),
        Line2D([], [], color="r", linestyle="--", linewidth=1.5, label="Schneider et al. 2020"),
        Line2D([], [], color="grey", linestyle="-", linewidth=1.5, label="NoSources Hot"),
    ]
    if include_nosources_cold:
        handles.append(Line2D([], [], color="k", linestyle="-", linewidth=1.5,
                              label="NoSources Cold"))
    return handles


def add_legends(ax, hot_cold_loc, runs_loc, ncol, handles, lines):
    first = ax.legend(handles=lines, loc=hot_cold_loc, fontsize=10, prop=LEGEND_FONT)
    ax.add_artist(first)
    ax.legend(handles=handles, loc=runs_loc, ncol=ncol, fontsize=10, prop=LEGEND_FONT)


def figure_velocity(fid, nodrag, nosrc):
    fig, ax = plt.subplots(figsize=(6.25, 5), dpi=300)
    r = fid["r"]
    lines = []
    for key, color, label in (("v_hot", "r", "Hot"), ("v_cold", "b", "Cold")):
        rr, median, p16, p84 = inner(r, *fid[key])
        lines += ax.plot(rr, median, lw=1.5, color=color, label=label)
        ax.fill_between(rr, p16, p84, alpha=0.1, color=color)
        ax.plot(*inner(r, nodrag[key][0]), lw=1.5, linestyle=":", color=color)
    ax.plot(*inner(r, nosrc["v_hot"][0]), lw=1.5, color="grey")
    ax.plot(*inner(r, nosrc["v_cold"][0]), lw=1.5, color="k")
    plot_schneider(ax, "vH", "vC")
    ax.axhline(0, linestyle="-", linewidth=0.7, alpha=0.7, color="k")

    add_legends(ax, "upper left", "upper right", 2, legend_runs(), lines)
    ax.set_xlabel("$R$   [kpc]")
    ax.set_ylabel("$v$  [km/s]")
    ax.set_ylim(-250, 2100)
    savefig(fig, "fig03_velocity.png")


def figure_density(fid, nodrag, nosrc):
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
    r = fid["r"]
    lines = []
    for key, color, label in (("n_hot", "r", "Hot"), ("n_cloud", "b", "Cold")):
        rr, median, p16, p84 = inner(r, *fid[key])
        lines += ax.plot(rr, median, lw=1.5, color=color, label=label)
        ax.fill_between(rr, p16, p84, alpha=0.1, color=color)
        ax.plot(*inner(r, nodrag[key][0]), lw=1.5, linestyle=":", color=color)
    # Without source terms no cold outflow exists: only the hot phase is shown.
    ax.plot(*inner(r, nosrc["n_hot"][0]), lw=1.5, color="grey")
    plot_schneider(ax, "nH", "nC")

    add_legends(ax, "upper right", "lower left", 1,
                legend_runs(include_nosources_cold=False), lines)
    ax.set_xlabel("$R$   [kpc]")
    ax.set_ylabel("$n$  [cm$^{-3}$]")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 100)
    savefig(fig, "fig04_density.png")


def figure_outflow_rate(fid, nodrag, nosrc):
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
    r = fid["r_mdot"]
    lines = []
    for key, color, label in (("mdot_hot", "r", "Hot"), ("mdot_cold", "b", "Cold")):
        rr, mean, std = inner(r, *fid[key])
        lines += ax.plot(rr, mean, lw=1.5, color=color, label=label)
        ax.fill_between(rr, mean - std, mean + std, alpha=0.1, color=color)
        ax.plot(*inner(r, nodrag[key][0]), lw=1.5, linestyle=":", color=color)
    ax.plot(*inner(r, nosrc["mdot_hot"][0]), lw=1.5, color="grey")
    ax.plot(*inner(r, nosrc["mdot_cold"][0]), lw=1.5, color="k")
    plot_schneider(ax, "MdotH", "MdotC")

    add_legends(ax, "upper left", "lower right", 2, legend_runs(), lines)
    ax.set_xlabel("$R$   [kpc]")
    ax.set_ylabel(r"$\dot{M}$   [M$_\odot$/yr]")
    ax.set_yscale("symlog", linthresh=0.3)
    ax.set_ylim(-0.17, 2)
    ax.yaxis.set_minor_locator(SymmetricalLogLocator(base=10, linthresh=0.3,
                                                     subs=[2, 3, 4, 5, 6, 7, 8, 9]))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="y", which="minor", length=2)
    savefig(fig, "fig05_outflow_rate.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sim-root", default=None)
    args = parser.parse_args()

    set_style(ticks_in=True)
    fid = time_averaged_run("fid", args.sim_root)
    nodrag = time_averaged_run("NoDrag_fmix0.03", args.sim_root)
    nosrc = time_averaged_run("NoSources", args.sim_root)

    figure_velocity(fid, nodrag, nosrc)
    figure_density(fid, nodrag, nosrc)
    figure_outflow_rate(fid, nodrag, nosrc)


if __name__ == "__main__":
    main()
