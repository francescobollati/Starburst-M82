"""
Figs. 6, 7, 8, 9, E.1, E.2 -- Effect of the mixing strength f_mix
(runs fid, fmix0.03, fmix0.1, fmix0.3, fmix1 and NoDrag_fmix0.03).

  Fig. 6   velocities and mass outflow rates of the hot and cold phases
  Fig. 7   source-term (and hot-phase pressure-gradient) accelerations
  Fig. 8   net hot-to-cold mass transfer rate, Mdot_net / dr
  Fig. 9   origin of the cold outflow: passive scalar P_cold,ic and Eq. (40)
  Fig. E.1 contributions to the hot-phase heating rate (fid vs fmix1)
  Fig. E.2 hot-phase temperature and total heating rate

Uses the pre-computed profiles in ``data/profiles/`` (see make_profiles.py).

Usage
-----
    python scripts/fig06_09_E1_E2_mixing.py
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from common import SERIF, Profiles, add_phase_titles, savefig, set_style

FMIX_RUNS = ["fid", "fmix0.03", "fmix0.1", "fmix0.3", "fmix1"]
FMIX_LABELS = ["0.003", "0.03", "0.1", "0.3", "1"]
COLORS = ["#253494", "#2c7fb8", "#41b6c4", "#7fcdbb", "#c7e9b4"]
LEGEND_FONT = {"family": SERIF, "size": 10}


def line_kwargs(i):
    """Colour per f_mix value; the fiducial run is drawn thicker."""
    kw = {"color": COLORS[i]}
    if i == 0:
        kw["linewidth"] = 3
    return kw


def fmix_legend(ax, loc, ncol=2, labels=None, handles=None, title=r"$f_{mix}$", **kw):
    return ax.legend(handles=handles, labels=labels, title=title, loc=loc, ncol=ncol,
                     fontsize=10, title_fontsize=10, **kw)


def fmt(x):
    return str(round(x, 2))


def total_heating(sg):
    """Net hot-phase heating rate: udot_mix + udot_drag - udot_cool."""
    return sg[5] + sg[6] - sg[11]


# ---------------------------------------------------------------------------
def fig6_velocities_and_rates(runs):
    fig, axs = plt.subplots(2, 2, figsize=(12, 10),
                            gridspec_kw={"wspace": 0.4, "hspace": 0.2}, dpi=300)
    for i, p in enumerate(runs):
        axs[0, 0].plot(p.velocities[0], p.velocities[1], label=FMIX_LABELS[i], **line_kwargs(i))
        axs[0, 1].plot(p.velocities[0], p.velocities[2], **line_kwargs(i))
        axs[1, 0].plot(p.mdots[0], p.mdots[1], **line_kwargs(i))
        axs[1, 1].plot(p.mdots[0], p.mdots[2], **line_kwargs(i))

    fmix_legend(axs[0, 0], "lower right")
    axs[0, 0].set_ylim(70, 1400)
    axs[0, 1].set_ylim(100, 800)
    for ax in axs[0]:
        ax.set_ylabel("$v$  [km/s]")
    for ax in axs[1]:
        ax.set_ylabel(r"$\dot{M}$   [$M_\odot$/yr]")
        ax.set_xlabel("$R$   [kpc]")
    add_phase_titles(axs[0, 0], axs[0, 1])
    savefig(fig, "fig06_velocity_outflow_rate.png")


def fig7_accelerations(runs, rmin=0.5):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5),
                            gridspec_kw={"wspace": 0.4, "hspace": 0.2}, dpi=300)
    hot_lines, cold_labels = [], []
    for i, p in enumerate(runs):
        a = p.accelerations[:, p.accelerations[0] > rmin]
        r, a_hydro, a_hot_drag, a_hot_mix, a_cold_drag, a_cold_mix = a
        hot_lines += axs[0].plot(r, a_hot_drag + a_hot_mix, label=FMIX_LABELS[i], **line_kwargs(i))
        axs[0].plot(r, a_hydro, linestyle="--", **line_kwargs(i))
        axs[1].plot(r, a_cold_drag + a_cold_mix, **line_kwargs(i))
        cold_labels.append(fmt(np.mean(np.abs(a_cold_mix) / a_cold_drag)))

    axs[0].add_artist(fmix_legend(axs[0], "lower right", handles=hot_lines))
    axs[0].legend(handles=[
        Line2D([], [], linestyle="solid", color="k", linewidth=1, label="$a_{drag}+a_{mix}$"),
        Line2D([], [], linestyle="--", color="k", linewidth=1, label="$a_{hydro}$"),
    ], loc="upper right", fontsize=10)
    fmix_legend(axs[1], "upper right", handles=axs[1].lines, labels=cold_labels,
                title=r"$\langle a_{c,mix}/a_{c,drag}\rangle$")
    for ax in axs:
        ax.set_ylabel("$a$   [cm/s$^2$]")
        ax.set_xlabel("$R$   [kpc]")
    add_phase_titles(axs[0], axs[1])
    savefig(fig, "fig07_accelerations.png")


def fig8_net_mass_transfer(runs, dr=0.3):
    fig, ax = plt.subplots(figsize=(6.3, 5), dpi=300)
    lines, ratios = [], []
    for i, p in enumerate(runs):
        r, mdot_loss, mdot_growth = p.mdots[0], p.mdots[3], p.mdots[4]
        lines += ax.plot(r, (mdot_growth - mdot_loss) / dr, label=FMIX_LABELS[i], **line_kwargs(i))
        ratios.append(fmt(np.mean(np.abs(mdot_growth - mdot_loss) / mdot_growth)))

    ax.add_artist(fmix_legend(ax, "lower left", handles=lines, prop=LEGEND_FONT))
    fmix_legend(ax, "lower right", handles=lines, labels=ratios, prop=LEGEND_FONT,
                title=r"$\langle |\dot{M}_{net}/\dot{M}_{growth} |\rangle$")
    ax.set_yscale("log")
    ax.set_ylabel(r"$\dot{M}_{net}/\Delta r$  [$M_\odot$/yr kpc$^{-1}$]")
    ax.set_xlabel("$R$   [kpc]")
    savefig(fig, "fig08_net_mass_transfer.png")


def fig9_cold_gas_origin(runs, nodrag):
    fig, ax = plt.subplots(1, 2, figsize=(12, 5),
                           gridspec_kw={"wspace": 0.4, "hspace": 0.2}, dpi=300)

    # Left: fraction of cold mass that was cold in the initial conditions.
    for i, p in enumerate(runs):
        ax[0].plot(p.scalars[0], p.scalars[4], label=FMIX_LABELS[i], **line_kwargs(i))
        if i == 1:  # NoDrag_fmix0.03 listed right after fmix0.03
            ax[0].plot(nodrag.scalars[0], nodrag.scalars[4], linestyle="--", color="k",
                       label="0.03  (NoDrag)")
    ax[0].add_artist(fmix_legend(ax[0], "upper right", ncol=1, prop=LEGEND_FONT))

    # Composition of the cold phase for fmix1: initially hot (ih) / cold (ic) / SN ejecta (isn).
    fmix1 = runs[-1].scalars
    ic = Line2D([], [], linestyle="solid", color=COLORS[4], linewidth=1, label=r"$\mathcal{P}_{cold,ic}$")
    isn, = ax[0].plot(fmix1[0], fmix1[7], linestyle=":", color=COLORS[4], label=r"$\mathcal{P}_{cold,isn}$")
    ih, = ax[0].plot(fmix1[0], fmix1[1], linestyle="--", color=COLORS[4], label=r"$\mathcal{P}_{cold,ih}$")
    ax[0].legend(handles=[ih, ic, isn], loc="lower right", fontsize=10)
    ax[0].set_ylim(0, 1)
    ax[0].set_ylabel(r"$\mathcal{P}_{cold,ic}$")
    ax[0].set_xlabel("$R$ [kpc]")

    # Right: outflow rate of initially cold gas, Eq. (40).
    for i, p in enumerate(runs):
        ratio = np.mean(p.mdots[5] / p.mdots[2])
        ax[1].plot(p.mdots[0], p.mdots[5], label=fmt(ratio), **line_kwargs(i))
    ax[1].legend(loc="lower right", prop=LEGEND_FONT, ncol=2, title_fontsize=10,
                 title=r"$\langle \dot{M}_{IC,cold} / \dot{M}_{cold} \rangle$")
    ax[1].set_ylabel(r"$\dot{M}_{IC,cold}$ [M$_\odot$/yr]")
    ax[1].set_xlabel("$R$ [kpc]")
    savefig(fig, "fig09_cold_gas_origin.png", bbox_inches="tight")


def figE1_heating_terms(runs, rmin=0.5):
    fig, ax = plt.subplots(figsize=(7.3, 5), dpi=300)
    styles = [("solid", 1, r"$\dot{u}_{tot}$"), ("-.", 2, r"$\dot{u}_{cool}$"),
              ("--", 2, r"$\dot{u}_{mix}$"), (":", 2, r"$\dot{u}_{drag}$")]
    handles = {}
    for i in (0, 4):  # fid and fmix1
        sg = runs[i].subgrid[:, runs[i].subgrid[0] > rmin]
        terms = [total_heating(sg), -sg[11], sg[5], sg[6]]
        for (ls, lw, label), y in zip(styles, terms):
            line, = ax.plot(sg[0], y, linestyle=ls, linewidth=lw, color=COLORS[i], label=label)
            handles.setdefault(i, []).append(line)

    ax.add_artist(ax.legend(handles=handles[0], loc="upper right", fontsize=10, ncol=2))
    fid_entry = Line2D([], [], linestyle="solid", color=COLORS[0], linewidth=1, label="$0.003$")
    handles[4][0].set_label("$1$")
    ax.legend(handles=[fid_entry, handles[4][0]], loc="lower right", title=r"$f_{mix}$",
              fontsize=10, title_fontsize=10)
    ax.set_yscale("symlog", linthresh=0.01)
    ax.set_ylabel(r"$\dot{u}_{tot}$   [(cm/s)$^2$s$^{-1}$]")
    ax.set_xlabel("$R$   [kpc]")
    savefig(fig, "figE1_heating_terms.png")


def figE2_temperature_heating(runs):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5),
                            gridspec_kw={"wspace": 0.4, "hspace": 0.2}, dpi=300)
    for i, p in enumerate(runs):
        r, temperature, mach = p.temperatures
        axs[0].plot(r, temperature, label=fmt(np.min(mach)), **line_kwargs(i))
        axs[1].plot(p.subgrid[0], total_heating(p.subgrid), label=FMIX_LABELS[i], **line_kwargs(i))

    fmix_legend(axs[0], "lower left", title=r"min$(\mathcal{M}_{hot})$")
    fmix_legend(axs[1], "lower left")
    axs[0].set_yscale("log")
    axs[0].set_ylabel("$T$   [K]")
    axs[1].set_yscale("symlog", linthresh=0.01)
    axs[1].set_ylabel(r"$\dot{u}_{tot}$   [(cm/s)$^2$s$^{-1}$]")
    for ax in axs:
        ax.set_xlabel("$R$   [kpc]")
    savefig(fig, "figE2_temperature_heating.png")


def main():
    set_style()
    runs = [Profiles(run) for run in FMIX_RUNS]
    nodrag = Profiles("NoDrag_fmix0.03")

    fig6_velocities_and_rates(runs)
    fig7_accelerations(runs)
    fig8_net_mass_transfer(runs)
    fig9_cold_gas_origin(runs, nodrag)
    figE1_heating_terms(runs)
    figE2_temperature_heating(runs)


if __name__ == "__main__":
    main()
