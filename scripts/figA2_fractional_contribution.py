"""
Fig. A.2 -- Fractional contribution of clouds with mass in (m_cl, m_cl + dm)
to the total drag force of a power-law cloud population, as a function of
m_crit (gamma = 0, left) and of gamma (m_crit = 1e2 Msun, right).

Usage
-----
    python scripts/figA2_fractional_contribution.py
"""

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from cloud_model import mass_to_radius, rinv_drag
from common import MSUN, savefig, set_style

M_MIN, M_MAX = 1e-2 * MSUN, 1e9 * MSUN
RHO_CL = 1e-25          # [g cm^-3]; cancels out in the fraction
N_GRID = 100

LEVELS = np.logspace(-10, 0, 1000)
NORM = mcolors.LogNorm(vmin=1e-10, vmax=1)
TICKS = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1]
FRACTION_LABEL = r"$\frac{\dot{p}(m_{cl})\,dm}{\int \dot{p}(m_{cl})\,dm}$"


def drag_fraction(m_crit, gamma, m_cl, dm):
    """Drag force of clouds in (m_cl, m_cl + dm) over the total drag force.

    Single-cloud drag scales as r_cl^2 (Eq. 8), weighted by the mass function
    (Eq. 17) and divided by the population integral (Eq. 27).
    """
    weight = np.where(m_cl < m_crit, (m_cl / m_crit)**(-gamma) / m_crit**2, 1.0 / m_cl**2)
    total = rinv_drag(m_crit, M_MIN, M_MAX, gamma, RHO_CL)
    return mass_to_radius(m_cl, RHO_CL)**2 * dm / total * weight


def panel(fig, ax, x, y, fraction, label_x):
    c = ax.contourf(x, y, fraction, levels=LEVELS, cmap="Blues", norm=NORM)
    fig.colorbar(c, ax=ax).set_ticks(TICKS)
    fig.text(label_x, 0.2, FRACTION_LABEL, fontsize=16, color="k", ha="right", va="bottom")
    ax.set_xscale("log")
    ax.set_xlabel(r"$m_{cl}$  [$M_\odot$]")


def main():
    set_style()
    m_edges = np.logspace(-2, 9, N_GRID) * MSUN
    m_cl, dm = m_edges[:-1], np.diff(m_edges)
    m_crits = np.logspace(-2, 9, N_GRID) * MSUN
    gammas = np.linspace(-1.99, 1.99, N_GRID)

    fig, ax = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={"wspace": 0.3}, dpi=300)

    mcl_grid, mcrit_grid = np.meshgrid(m_cl, m_crits)
    fraction = drag_fraction(mcrit_grid, 0.0, mcl_grid, dm)
    panel(fig, ax[0], mcl_grid / MSUN, mcrit_grid / MSUN, fraction, 0.35)
    ax[0].set_yscale("log")
    ax[0].set_ylabel(r"$m_{crit}$  [$M_\odot$]")

    mcl_grid, gamma_grid = np.meshgrid(m_cl, gammas)
    fraction = drag_fraction(1e2 * MSUN, gamma_grid, mcl_grid, dm)
    panel(fig, ax[1], mcl_grid / MSUN, gamma_grid, fraction, 0.8)
    ax[1].set_ylabel(r"$\gamma$")

    savefig(fig, "figA2_fractional_contribution.png")


if __name__ == "__main__":
    main()
