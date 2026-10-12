"""
Analytic coefficients of the piecewise power-law cloud mass distribution
(Sect. 2.2 of the paper), dN/dm = N * (m/m_crit)^-gamma / m_crit^2 for
m < m_crit and N * m^-2 above, defined on (m_min, m_max).

All functions accept cgs masses [g] and densities [g cm^-3] and broadcast
over numpy arrays. ``gamma`` is the low-mass slope.

``clip=True`` replaces non-positive coefficients by a tiny positive number,
which is how the coefficients are regularised in the simulation code
(Sect. 2.2) and is used to draw them on a log scale (Fig. A.1).
"""

import numpy as np

_TINY = 1e-60


def mass_to_radius(mass, rho):
    return (mass / (4.0 / 3.0 * np.pi * rho))**(1.0 / 3.0)


def radius_to_mass(radius, rho):
    return 4.0 / 3.0 * np.pi * rho * radius**3


def _clip(values, clip):
    return np.where(values > 0, values, _TINY) if clip else values


def normalisation(m_crit, m_min, m_max, gamma):
    """N / M_cold, Eq. (18)."""
    x = m_min / m_crit
    return 1.0 / ((1.0 - x**(2.0 - gamma)) / (2.0 - gamma) + np.log(m_max / m_crit))


def rinv_loss(m_crit, m_min, m_max, gamma, rho_cl, clip=False):
    """r_loss^-1, Eq. (21) [cm^-1]."""
    x = m_min / m_crit
    r_crit, r_min, r_max = (mass_to_radius(m, rho_cl) for m in (m_crit, m_min, m_max))
    value = (9.0 * (2.0 - gamma) / (5.0 - 3.0 * gamma) / r_crit
             - 3.0 / (5.0 - 3.0 * gamma) * x**(2.0 - gamma) / r_min
             - 3.0 / r_max)
    return _clip(value, clip)


def rinv_growth(m_crit, m_min, m_max, gamma, rho_cl, clip=False):
    """r_growth^-1, Eq. (23) [cm^-1]."""
    x = m_min / m_crit
    r_crit, r_min, r_max = (mass_to_radius(m, rho_cl) for m in (m_crit, m_min, m_max))
    xi_min, xi_max = r_min / r_crit, r_max / r_crit
    value = (6.0 / (11.0 - 6.0 * gamma) * (1.0 / r_crit - x**(2.0 - gamma) * np.sqrt(xi_min) / r_min)
             - 4.0 * (xi_max**0.25 / r_max - 1.0 / r_crit))
    return _clip(value, clip)


def rinv_drag(m_crit, m_min, m_max, gamma, rho_cl, clip=False):
    """r_drag^-1 / rho, Eq. (27) without the density prefactor [cm^2 g^-1]."""
    x = m_min / m_crit
    r_crit, r_min, r_max = (mass_to_radius(m, rho_cl) for m in (m_crit, m_min, m_max))
    value = (3.0 / (5.0 - 3.0 * gamma) * (r_crit**2 / m_crit - x**(2.0 - gamma) * r_min**2 / m_min)
             - 3.0 * (r_max**2 / m_max - r_crit**2 / m_crit))
    return _clip(value, clip)
