"""
Compute the time-averaged radial profiles of one simulation and store them as
``.npy`` files. These files are the input of the profile-based figure scripts
(Figs. 6-9, D.1, E.1, E.2, F.1). The profiles used for the paper are already
provided in ``data/profiles/``; run this script only to regenerate them from
the snapshots.

Usage
-----
    python scripts/make_profiles.py fid
    python scripts/make_profiles.py fmix0.1 --sim-root /path/to/Simulations

The layout of each output array is documented in ``common.Profiles``.
"""

import argparse
import os

import numpy as np

from common import (
    DR, GAMMA, K_B, M_P, MU, MSUN, T_COLD, UNIT_DENS, UNIT_FORCE, UNIT_LEN,
    UNIT_MASS, UNIT_PRES, UNIT_TIME, UNIT_VEL, YR, PROFILE_DIR, RUN_DIRS,
    Profiles, SliceGrid, averaging_snapshots, outflow_rate, radial_profile, run_output_dir,
    shell_sum,
)

N_PIXELS = 128       # resolution of the averaging slice
RMAX_PROFILE = 5.5   # outer radius of the profiles [kpc]
RMAX_RATES = 8.75    # outer radius used when binning the outflow rates [kpc]


def cooling_rate(temperature, density_cgs):
    """Radiative cooling rate du/dt [erg s^-1 g^-1] of the hot phase.

    Piecewise parabolic fit to a solar-metallicity cooling curve
    (Schneider & Robertson 2018), as used in the simulations.
    Note: for log T < 4 the exponent is 0 (as in the original analysis);
    the hot phase never reaches such temperatures in the outflow cone.
    """
    logt = np.log10(temperature)
    exponent = np.zeros_like(logt)
    exponent = np.where((logt >= 4.0) & (logt < 5.9), -1.3 * (logt - 5.25)**2 - 21.25, exponent)
    exponent = np.where((logt >= 5.9) & (logt < 7.4), 0.7 * (logt - 7.1)**2 - 22.8, exponent)
    exponent = np.where(logt >= 7.4, 0.45 * logt - 26.065, exponent)
    n_h = 0.76 * density_cgs / MU / M_P
    return n_h * n_h * 10.0**exponent / density_cgs


def cell_quantities(snap):
    """All per-cell quantities that are averaged on the slice."""
    vel_hot, vel_cold = snap["Velocities"], snap["plVelocities"]
    sources = snap["SourcesOutputCloud"]
    mdot_growth = sources[:, 7]          # code mass / code time
    mdot_loss = sources[:, 8]
    pdot_drag_code = sources[:, 12]

    # Thermodynamics of the hot phase
    u_hot = snap["InternalEnergy"]
    u_cold = T_COLD * K_B / MU / M_P / (GAMMA - 1.0) / UNIT_VEL**2
    sound_speed = np.sqrt(GAMMA * snap.pressure * UNIT_PRES / (snap.density * UNIT_DENS)) / UNIT_VEL
    mach = np.linalg.norm(vel_hot, axis=1) / sound_speed

    grad_p = snap["PressureGradient"]
    grad_p_norm = np.linalg.norm(grad_p, axis=1)
    p_scale = np.zeros_like(snap.density)
    p_scale[grad_p_norm > 0] = snap.pressure[grad_p_norm > 0] / grad_p_norm[grad_p_norm > 0]
    pressure_force = snap.radial(grad_p) * snap.volume * UNIT_FORCE

    # Energy source terms of the hot phase, Eqs. (25) and (28) [cm^2 s^-3]
    vrel_3d = np.linalg.norm(vel_hot - vel_cold, axis=1)
    udot_mix = mdot_loss / snap.cold_mass * (u_cold - u_hot + 0.5 * vrel_3d**2) * UNIT_VEL**2 / UNIT_TIME
    udot_drag = pdot_drag_code * vrel_3d / snap.mass * UNIT_VEL**2 / UNIT_TIME
    udot_ratio = np.where(udot_drag > 0, udot_mix / udot_drag, 0.0)

    # Source-term accelerations, Eqs. (36)-(38) [cm s^-2]
    pdot_mix = sources[:, 11] * UNIT_FORCE
    pdot_drag = pdot_drag_code * UNIT_FORCE
    vrel_rad = snap.vrad - snap.cold_vrad
    hot_mass_cgs, cold_mass_cgs = snap.mass * UNIT_MASS, snap.cold_mass * UNIT_MASS

    return {
        "vrad_hot": snap.vrad,
        "vrad_cold": snap.cold_vrad,
        "n_hot": snap.density * UNIT_DENS / MU / M_P,
        "n_cloud": snap.cloud_density * UNIT_DENS / MU / M_P,
        "temperature": snap.temperature,
        "mach": mach,
        "p_scale": p_scale,
        "cell_size": snap.volume**(1.0 / 3.0),
        "r_cl": snap.cloud_radius,
        "m_cl": snap.cloud_mass * UNIT_MASS / MSUN,
        "pdot_mix": pdot_mix,
        "pdot_drag": pdot_drag,
        "udot_mix": udot_mix,
        "udot_drag": udot_drag,
        "t_cool": sources[:, 4],
        "udot_ratio": udot_ratio,
        "udot_cool": cooling_rate(snap.temperature, snap.density * UNIT_DENS),
        "pressure_force": pressure_force,
        "a_hydro": -pressure_force / hot_mass_cgs,
        "a_hot_drag": -pdot_drag / hot_mass_cgs,
        "a_hot_mix": -(mdot_loss / snap.mass / UNIT_TIME) * vrel_rad * UNIT_VEL,
        "a_cold_drag": pdot_drag / cold_mass_cgs,
        "a_cold_mix": (mdot_growth / snap.cold_mass / UNIT_TIME) * vrel_rad * UNIT_VEL,
        **{f"P_cold_{i}": snap["plPassiveScalars"][:, k] for k, i in enumerate(("ih", "ic", "isn"))},
        **{f"P_hot_{i}": snap["PassiveScalars"][:, k] for k, i in enumerate(("ih", "ic", "isn"))},
    }


def shell_rates(snap):
    """Mass fluxes through radial shells inside the outflow cone [Msun/yr]."""
    sources = snap["SourcesOutputCloud"]
    to_msun_yr = UNIT_MASS / UNIT_TIME * YR / MSUN
    radii, mdot_hot = outflow_rate(snap, snap.mass, snap.vrad, RMAX_RATES)
    _, mdot_cold = outflow_rate(snap, snap.cold_mass, snap.cold_vrad, RMAX_RATES)
    _, mdot_loss = shell_sum(snap, sources[:, 8], RMAX_RATES)
    _, mdot_growth = shell_sum(snap, sources[:, 7], RMAX_RATES)

    # Outflow rate of initially cold gas, Eq. (40)
    flux_hot_ic = snap.mass * snap["PassiveScalars"][:, 1] * snap.vrad
    flux_cold_ic = snap.cold_mass * snap["plPassiveScalars"][:, 1] * snap.cold_vrad
    _, sum_hot_ic = shell_sum(snap, flux_hot_ic, RMAX_RATES)
    _, sum_cold_ic = shell_sum(snap, flux_cold_ic, RMAX_RATES)
    mdot_ic_cold = (sum_hot_ic + sum_cold_ic) / DR * UNIT_MASS * UNIT_VEL / UNIT_LEN * YR / MSUN

    return radii, np.array([mdot_hot, mdot_cold, mdot_loss * to_msun_yr,
                            mdot_growth * to_msun_yr, mdot_ic_cold,
                            sum_cold_ic / sum_hot_ic])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run", choices=sorted(RUN_DIRS), help="run name as in Table 2")
    parser.add_argument("--sim-root", default=None,
                        help="folder containing the simulation directories "
                             "(default: $M82_SIM_DIR or ./Simulations)")
    parser.add_argument("--out", default=PROFILE_DIR, help="output root folder")
    args = parser.parse_args()

    grid, maps, rates = None, {}, []
    n_snap = 0
    for snap in averaging_snapshots(run_output_dir(args.run, args.sim_root)):
        print(f"{os.path.basename(snap.filename)}  t = {snap.time_myr:.2f} Myr")
        if grid is None:
            grid = SliceGrid(snap.boxsize, N_PIXELS)
        idx = grid.nearest_cells(snap)
        for key, values in cell_quantities(snap).items():
            maps[key] = maps.get(key, 0.0) + grid.sample(values, idx)
        rate_radii, snap_rates = shell_rates(snap)
        rates.append(snap_rates)
        n_snap += 1

    def profile(key, *stats):
        radii, median, p16, p84 = radial_profile(grid, maps[key] / n_snap, RMAX_PROFILE)
        chosen = {"median": median, "p16": p16, "p84": p84}
        return [chosen[s] for s in (stats or ("median",))]

    radii = radial_profile(grid, maps["vrad_hot"], RMAX_PROFILE)[0]
    mean_rates = np.mean(rates, axis=0)[:, rate_radii < RMAX_PROFILE]

    outputs = {
        "velocities": [radii, *profile("vrad_hot"), *profile("vrad_cold")],
        "densities": [radii, *profile("n_hot"), *profile("n_cloud")],
        "temperatures": [radii, *profile("temperature"), *profile("mach")],
        "resolution": [radii, *profile("p_scale"), *profile("cell_size")],
        "mdots": [rate_radii[rate_radii < RMAX_PROFILE], *mean_rates],
        "subgrid": [radii] + [p for key in ("r_cl", "m_cl", "pdot_mix", "pdot_drag",
                                            "udot_mix", "udot_drag", "t_cool")
                              for p in profile(key)]
                   + profile("udot_ratio", "median", "p16", "p84")
                   + profile("udot_cool") + profile("pressure_force"),
        "accelerations": [radii] + [p for key in ("a_hydro", "a_hot_drag", "a_hot_mix",
                                                  "a_cold_drag", "a_cold_mix")
                                    for p in profile(key)],
        "scalars": [radii] + [p for key in ("P_cold_ih", "P_cold_ic", "P_cold_isn",
                                            "P_hot_ih", "P_hot_ic", "P_hot_isn")
                              for p in profile(key, "median", "p16", "p84")],
    }

    out_dir = os.path.join(args.out, args.run)
    os.makedirs(out_dir, exist_ok=True)
    for name, rows in outputs.items():
        np.save(os.path.join(out_dir, Profiles.FILES[name]), np.array(rows))
    print(f"profiles written to {out_dir}")


if __name__ == "__main__":
    main()
