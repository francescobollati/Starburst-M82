"""
Fig. D.1 -- Resolution test: radial velocity profiles of the hot (left) and
cold (right) phases for the Vc* runs, which differ from fid only in the target
cell volume V_crit of the outflow region.

Uses the pre-computed profiles in ``data/profiles/`` (see make_profiles.py).

Usage
-----
    python scripts/figD1_resolution.py
"""

import matplotlib.pyplot as plt

from common import SERIF, Profiles, add_phase_titles, savefig, set_style

RUNS = ["Vc0.001", "Vc0.01", "Vc0.03", "fid", "Vc0.3", "Vc1"]
LABELS = ["$0.001$", "$0.01$", "$0.03$", "$0.1$ (fid)", "$0.3$", "$1$"]
COLORS = ["#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177"]


def main():
    set_style()
    fig, axs = plt.subplots(1, 2, figsize=(12, 5),
                            gridspec_kw={"wspace": 0.4, "hspace": 0.2}, dpi=300)
    for run, label, color in zip(RUNS, LABELS, COLORS):
        r, v_hot, v_cold = Profiles(run).velocities
        width = 3 if run == "fid" else None
        axs[0].plot(r, v_hot, color=color, linewidth=width, label=label)
        axs[1].plot(r, v_cold, color=color, linewidth=width)

    axs[0].legend(prop={"family": SERIF, "size": 10}, title=r"$V_{crit}$  [kpc$^3$]",
                  title_fontproperties={"family": SERIF, "size": 10})
    axs[0].set_ylim(250, 1600)
    for ax in axs:
        ax.set_xlabel("$R$   [kpc]")
        ax.set_ylabel("$v$  [km/s]")
    add_phase_titles(axs[0], axs[1])
    savefig(fig, "figD1_resolution.png")


if __name__ == "__main__":
    main()
