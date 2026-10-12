"""
Shared constants, unit conversions, snapshot I/O and analysis helpers used by
all figure scripts of Bollati & Weinberger (2026, A&A 710, A88).

Conventions
-----------
* AREPO code units: length = 1 kpc, velocity = 1 km/s, mass = 1e10 Msun.
* The galactic disc lies in the plane z = Boxsize/2, so the outflow
  propagates along +-z.
* "Hot" quantities refer to the gas phase evolved by the hydro solver;
  "cold" quantities (prefix ``pl`` in the snapshots) refer to the
  pressureless cold phase.
* Time-averaged quantities are computed over snapshots 119-128, i.e. the
  interval T = (27, 30) Myr used throughout Sect. 4 of the paper.
* Radial profiles are restricted to a bi-cone of semi-opening angle 30 deg
  around the z axis and use radial bins of 0.3 kpc.
"""

import json
import os
from functools import cached_property

import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# Physical constants (cgs). Values kept identical to those used for the paper.
# ---------------------------------------------------------------------------
GAMMA = 5.0 / 3.0
K_B = 1.3806488e-16
M_P = 1.67373522381e-24
MU = 0.6                      # mean molecular weight
YR = 3.1536e7
MYR = 3.1536e13
PC = 3.086e18
KPC = 1.0e3 * PC
MSUN = 2.0e33
T_COLD = 1.0e4                # temperature of the cold phase [K]

# ---------------------------------------------------------------------------
# AREPO code units -> cgs
# ---------------------------------------------------------------------------
UNIT_VEL = 1.0e5
UNIT_LEN = 3.085678e21
UNIT_MASS = 1.989e43
UNIT_TIME = UNIT_LEN / UNIT_VEL
UNIT_DENS = UNIT_MASS / UNIT_LEN**3
UNIT_PRES = UNIT_VEL**2 * UNIT_DENS
UNIT_FORCE = UNIT_PRES * UNIT_LEN**2

# ---------------------------------------------------------------------------
# Analysis choices (Sect. 4 of the paper)
# ---------------------------------------------------------------------------
SNAP_FIRST, SNAP_LAST = 119, 128        # time-averaging window, T = (27, 30) Myr
CONE_HALF_ANGLE = np.pi / 6.0           # 30 deg bi-cone around the outflow axis
DR = 0.3                                # radial bin width [kpc]

# ---------------------------------------------------------------------------
# Paths and run names
# ---------------------------------------------------------------------------
REPO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(REPO_DIR, "data")
PROFILE_DIR = os.path.join(DATA_DIR, "profiles")
FIGURE_DIR = os.path.join(REPO_DIR, "figures")

# Run names used in the paper (Table 2) -> directory names of the simulations.
RUN_DIRS = {
    "fid": "fid",
    "NoSources": "15J_NoSources",
    "NoDrag_fmix0.03": "fc0.03_nd",
    "fmix0.03": "fc0.03",
    "fmix0.1": "fc0.1",
    "fmix0.3": "fc0.3",
    "fmix1": "fc1",
    "Vc0.001": "high_res_rainer",
    "Vc0.01": "HR1",
    "Vc0.03": "HR_Vc0.03",
    "Vc0.3": "HR_Vc0.3",
    "Vc1": "HR_Vc1",
    "fcold0.1": "fcld_l",
    "fcold0.9": "fcld_h",
    "rcl0.01": "rcl0.01",
    "rcl0.6": "rcl0.6",
}

SERIF = "Times New Roman"   # font used for legends and panel titles


def default_sim_root():
    """Root directory containing the simulation folders (see ``RUN_DIRS``)."""
    return os.environ.get("M82_SIM_DIR", "Simulations")


def run_output_dir(run, sim_root=None):
    """Directory holding the ``snap_XXX.hdf5`` files of a run."""
    return os.path.join(sim_root or default_sim_root(), RUN_DIRS[run], "output")


def savefig(fig, name, **kwargs):
    """Save a figure into ``figures/`` and report where it went."""
    os.makedirs(FIGURE_DIR, exist_ok=True)
    path = os.path.join(FIGURE_DIR, name)
    fig.savefig(path, **kwargs)
    print(f"saved {os.path.relpath(path, REPO_DIR)}")


# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------
def set_style(ticks_in=False):
    """Matplotlib style of the paper figures.

    ``ticks_in=True`` reproduces the style of Figs. 2-5 (inward ticks on all
    four sides with minor ticks); the other figures use matplotlib defaults.
    """
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
    plt.rcParams.update({
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
    })
    if ticks_in:
        plt.rcParams.update({
            "xtick.direction": "in", "ytick.direction": "in",
            "xtick.top": True, "ytick.right": True,
            "xtick.minor.visible": True, "ytick.minor.visible": True,
            "lines.dash_capstyle": "round", "lines.solid_capstyle": "round",
            "legend.handletextpad": 0.4, "legend.handlelength": 2,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6, "ytick.major.width": 0.6,
            "xtick.minor.width": 0.45, "ytick.minor.width": 0.45,
            "xtick.major.size": 2.75, "ytick.major.size": 2.75,
            "xtick.minor.size": 1.75, "ytick.minor.size": 1.75,
            "figure.dpi": 200,
        })


def add_phase_titles(ax_hot, ax_cold):
    """'Hot phase' / 'Cold phase' headers above a pair of panels."""
    for ax, title in ((ax_hot, "Hot phase"), (ax_cold, "Cold phase")):
        ax.text(0.5, 1.05, title, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=12, fontweight="bold", fontname=SERIF)


# ---------------------------------------------------------------------------
# Snapshot reading
# ---------------------------------------------------------------------------
class Snapshot:
    """Lazy access to the gas cells of one AREPO multifluid snapshot.

    Raw fields are in code units; derived quantities are documented in
    their docstrings. Positions are measured from the box centre.
    """

    def __init__(self, filename):
        self.filename = filename
        with h5py.File(filename, "r") as f:
            self.boxsize = float(f["Header"].attrs["BoxSize"])
            self.time_myr = float(f["Header"].attrs["Time"]) * UNIT_TIME / MYR
            self._raw = {key: np.array(f["PartType0"][key], dtype=np.float64)
                         for key in f["PartType0"].keys()
                         if key not in ("ParticleIDs",)}

    def __getitem__(self, key):
        return self._raw[key]

    # --- geometry ---------------------------------------------------------
    @cached_property
    def pos(self):
        """Cell positions relative to the box centre [kpc]."""
        return self["Coordinates"] - 0.5 * self.boxsize

    @cached_property
    def radius(self):
        """Spherical radius [kpc]."""
        return np.linalg.norm(self.pos, axis=1)

    @cached_property
    def theta(self):
        """Angle from the outflow (z) axis, folded to [0, pi/2]."""
        return np.arccos(np.abs(self.pos[:, 2]) / self.radius)

    def radial(self, vec):
        """Radial component of a vector field (spherical radius)."""
        return np.sum(vec * self.pos, axis=1) / self.radius

    # --- hot phase --------------------------------------------------------
    @property
    def mass(self):
        return self["Masses"]

    @property
    def density(self):
        return self["Density"]

    @property
    def volume(self):
        return self["Volume"]

    @cached_property
    def vrad(self):
        """Radial velocity of the hot phase [km/s]."""
        return self.radial(self["Velocities"])

    @cached_property
    def pressure(self):
        """Hot-phase pressure [code units]."""
        return (GAMMA - 1.0) * self.density * self["InternalEnergy"]

    @cached_property
    def temperature(self):
        """Hot-phase temperature [K]."""
        return (self.pressure * UNIT_PRES / (self.density * UNIT_DENS)
                * MU * M_P / K_B)

    # --- cold (pressureless) phase ----------------------------------------
    @property
    def cold_density(self):
        """Cold mass per cell volume [code units]."""
        return self["plDensity"]

    @cached_property
    def cold_mass(self):
        return self.cold_density * self.volume

    @cached_property
    def cold_vrad(self):
        """Radial velocity of the cold phase [km/s]."""
        return self.radial(self["plVelocities"])

    # --- sub-grid cloud model output (SourcesOutputCloud columns) ---------
    @property
    def cloud_radius(self):
        """r_cl (delta model) or r_crit (power-law model) [kpc]."""
        return self["SourcesOutputCloud"][:, 0]

    @property
    def cloud_mass(self):
        """m_cl (delta model) or m_crit (power-law model) [code mass]."""
        return self["SourcesOutputCloud"][:, 1]

    @property
    def chi(self):
        """Cloud-to-hot density contrast."""
        return self["SourcesOutputCloud"][:, 2]

    @property
    def cloud_density(self):
        """Density of the individual sub-grid clouds, Eq. (11) [code units]."""
        return self.chi * self.density


def snapshot_path(output_dir, i_snap):
    return os.path.join(output_dir, "snap_%03d.hdf5" % i_snap)


def averaging_snapshots(output_dir, first=SNAP_FIRST, last=SNAP_LAST):
    """Yield the snapshots of the time-averaging window."""
    for i_snap in range(first, last + 1):
        yield Snapshot(snapshot_path(output_dir, i_snap))


# ---------------------------------------------------------------------------
# Edge-on slice through the box centre
# ---------------------------------------------------------------------------
class SliceGrid:
    """Regular ``n x n`` grid on the x-z plane through the box centre.

    Each pixel takes the value of the nearest gas cell. Arrays returned by
    :meth:`sample` have shape ``(n, n)``, with the first index running along
    z (vertical) and the second along x (horizontal).
    """

    def __init__(self, boxsize, n):
        self.n = n
        self.boxsize = boxsize
        self.edges = np.linspace(0.0, boxsize, n + 1)
        centres = 0.5 * (self.edges[1:] + self.edges[:-1])
        xx, zz = np.meshgrid(centres, centres)
        self.points = np.column_stack(
            [xx.ravel(), np.full(n * n, 0.5 * boxsize), zz.ravel()])
        dx, dz = xx - 0.5 * boxsize, zz - 0.5 * boxsize
        self.radius = np.hypot(dx, dz)
        self.theta = np.arccos(np.abs(dz) / self.radius)
        self.in_cone = self.theta <= CONE_HALF_ANGLE

    def nearest_cells(self, snap):
        """Index of the gas cell closest to each pixel."""
        _, idx = cKDTree(snap["Coordinates"]).query(self.points, k=1)
        return idx

    def sample(self, values, idx):
        return values[idx].reshape(self.n, self.n)


def radial_profile(grid, image, rmax, dr=DR):
    """Median and 16th/84th percentiles of ``image`` in radial shells.

    Only pixels inside the outflow cone are used. Returns
    ``(radii, median, p16, p84)`` with ``radii`` the bin centres [kpc].
    """
    nbins = int(rmax / dr)
    radii = (np.arange(nbins) + 0.5) * dr
    median, p16, p84 = (np.empty(nbins) for _ in range(3))
    for i in range(nbins):
        sel = grid.in_cone & (grid.radius >= i * dr) & (grid.radius < (i + 1) * dr)
        median[i], p16[i], p84[i] = np.percentile(image[sel], [50, 16, 84])
    return radii, median, p16, p84


def shell_sum(snap, values, rmax, dr=DR):
    """Sum of a per-cell quantity in radial shells inside the outflow cone.

    Returns ``(radii, sums)``. Used for the mass outflow rate (Sect. 4.1):
    ``Mdot(r) = sum_i m_i v_i / dr``.
    """
    nbins = int(rmax / dr)
    radii = (np.arange(nbins) + 0.5) * dr
    cone = snap.theta <= CONE_HALF_ANGLE
    sums = np.empty(nbins)
    for i in range(nbins):
        sel = cone & (snap.radius >= i * dr) & (snap.radius < (i + 1) * dr)
        sums[i] = np.sum(values[sel])
    return radii, sums


def outflow_rate(snap, mass, vrad, rmax, dr=DR):
    """Mass outflow rate in radial shells inside the cone [Msun/yr]."""
    radii, flux = shell_sum(snap, mass * vrad, rmax, dr)
    return radii, flux / dr * UNIT_MASS * UNIT_VEL / UNIT_LEN * YR / MSUN


# ---------------------------------------------------------------------------
# Pre-computed radial profiles (written by make_profiles.py)
# ---------------------------------------------------------------------------
class Profiles:
    """Time-averaged radial profiles of one run, as saved by ``make_profiles.py``.

    Every array has the radius [kpc] in row 0. Row meanings:

    velocities    : vrad_hot, vrad_cold                                  [km/s]
    densities     : n_hot, n_cloud (sub-grid cloud density, Eq. 11)     [cm^-3]
    temperatures  : T_hot [K], Mach_hot
    resolution    : pressure scale height, cell size                     [kpc]
    mdots         : Mdot_hot, Mdot_cold, Mdot_loss, Mdot_growth,
                    Mdot_IC,cold (Eq. 40), ratio of cold/hot IC-cold flux [Msun/yr]
    subgrid       : r_crit [kpc], m_crit [Msun], pdot_mix, pdot_drag [dyn],
                    udot_mix, udot_drag [cm^2 s^-3], t_cool [code],
                    udot_mix/udot_drag (median, p16, p84),
                    udot_cool [cm^2 s^-3], pressure force [dyn]
    accelerations : a_hydro, a_drag, a_mix (hot); a_drag, a_mix (cold)   [cm s^-2]
    scalars       : median, p16, p84 of P_cold,ih, P_cold,ic, P_cold,isn,
                    P_hot,ih, P_hot,ic, P_hot,isn
    """

    FILES = {
        "velocities": "velocities.npy",
        "densities": "densities.npy",
        "temperatures": "Temperatures.npy",
        "resolution": "Resolution.npy",
        "mdots": "Mdots.npy",
        "subgrid": "subgrid.npy",
        "accelerations": "accelerations.npy",
        "scalars": "scalars.npy",
    }

    def __init__(self, run, profile_dir=PROFILE_DIR):
        self.run = run
        self.directory = os.path.join(profile_dir, run)

    def __getattr__(self, name):
        if name not in self.FILES:
            raise AttributeError(name)
        data = np.load(os.path.join(self.directory, self.FILES[name]))
        setattr(self, name, data)
        return data


# ---------------------------------------------------------------------------
# Digitised data from Schneider et al. (2020)
# ---------------------------------------------------------------------------
def read_digitized(name):
    """Curve digitised (WebPlotDigitizer) from Schneider et al. (2020).

    Returns ``(r, y)`` sorted by radius.
    """
    with open(os.path.join(DATA_DIR, "schneider2020", name + ".json")) as f:
        points = json.load(f)
    r = np.array([float(p["x"]) for p in points])
    y = np.array([float(p["y"]) for p in points])
    order = np.argsort(r)
    return r[order], y[order]
