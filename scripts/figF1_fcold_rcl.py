"""
Fig. F.1 -- Dependence on the initial cold gas fraction f_cold (fcold* runs)
and comparison with the delta-distribution cloud model (rcl* runs).

Panels: cold mass outflow rate, cold-phase velocity, and characteristic cloud
radius (r_cl for the rcl* runs, r_crit for all the others).

Uses the pre-computed profiles in ``data/profiles/`` (see make_profiles.py).

Usage
-----
    python scripts/figF1_fcold_rcl.py
"""

import matplotlib.pyplot as plt

from common import SERIF, Profiles, savefig, set_style

# run name, legend label, colour, line style
RUNS = [
    ("fid", "fid", "k", "-"),
    ("fcold0.1", "$f_{cold}$ = 0.1", "mediumturquoise", "-"),
    ("fcold0.9", "$f_{cold}$ = 0.9", "mediumturquoise", "--"),
    ("rcl0.01", "$r_{cl}$ = 0.01 kpc", "sandybrown", "-"),
    ("rcl0.6", "$r_{cl}$ = 0.6 kpc", "sandybrown", "--"),
]
RMIN_RADIUS = 0.5   # inner radius of the cloud-radius panel [kpc]


def main():
    set_style()
    fig, axs = plt.subplots(1, 3, figsize=(16, 5),
                            gridspec_kw={"wspace": 0.4, "hspace": 0.2}, dpi=300)
    for run, label, color, ls in RUNS:
        p = Profiles(run)
        kw = dict(color=color, linestyle=ls, label=label,
                  linewidth=3 if run == "fid" else None)
        axs[0].plot(p.mdots[0], p.mdots[2], **kw)
        axs[1].plot(p.velocities[0], p.velocities[2], **kw)
        outer = p.subgrid[0] > RMIN_RADIUS
        axs[2].plot(p.subgrid[0, outer], p.subgrid[1, outer], **kw)

    axs[0].set_ylabel(r"$\dot{M}_{cold}$  [M$_\odot$/yr]")
    axs[1].legend(prop={"family": SERIF, "size": 10}, loc="upper left")
    axs[1].set_ylim(-100, 800)
    axs[1].set_ylabel("$v_{cold}$  [km/s]")
    axs[2].set_yscale("log")
    axs[2].set_ylim(1e-4, 1)
    axs[2].set_ylabel("$r$  [kpc]")
    for ax in axs:
        ax.set_xlabel("$R$   [kpc]")
    savefig(fig, "figF1_fcold_rcl.png")


if __name__ == "__main__":
    main()
