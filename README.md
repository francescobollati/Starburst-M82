# Starburst-M82

Analysis and plotting scripts for

> **Modeling multiphase galactic outflows: A multifluid moving-mesh approach**
> F. Bollati & R. Weinberger, *Astronomy & Astrophysics* **710**, A88 (2026)
> [doi:10.1051/0004-6361/202556959](https://doi.org/10.1051/0004-6361/202556959)

The paper presents a sub-grid model for multiphase galactic outflows in the
multifluid version of the moving-mesh code AREPO. The cold phase (T = 10⁴ K)
is a pressureless fluid made of unresolved clouds, and it exchanges mass,
momentum and energy with the hot phase through drag and mixing source terms
calibrated on cloud-crushing simulations. We apply the model to an
M82-like starburst galaxy and compare the results with the high-resolution
simulation of Schneider et al. (2020).

This repository contains the code used to produce every figure of the paper
from the simulation outputs, together with the processed data needed to
remake most of them without the snapshots.

## Figures and scripts

| Figure | Content | Script | Input |
|---|---|---|---|
| 1 | Schematic of the model | – (made by hand) | – |
| 2 | Density slices, NoSources vs fid | [`fig02_density_slices.py`](scripts/fig02_density_slices.py) | snapshots |
| 3, 4, 5 | Velocity, density and mass outflow rate vs Schneider et al. (2020) | [`fig03_04_05_schneider_comparison.py`](scripts/fig03_04_05_schneider_comparison.py) | snapshots + `data/schneider2020` |
| 6 | Velocities and outflow rates vs *f*<sub>mix</sub> | [`fig06_09_E1_E2_mixing.py`](scripts/fig06_09_E1_E2_mixing.py) | `data/profiles` |
| 7 | Source-term accelerations | ″ | `data/profiles` |
| 8 | Net hot→cold mass transfer | ″ | `data/profiles` |
| 9 | Origin of the cold outflow (passive scalars) | ″ | `data/profiles` |
| A.1 | r<sub>drag</sub>⁻¹(m<sub>crit</sub>) for varying γ, m<sub>min</sub>, m<sub>max</sub> | [`figA1_rinv_coefficients.py`](scripts/figA1_rinv_coefficients.py) | analytic + `data/fid_mcrit_histogram.npz` |
| A.2 | Fractional contribution of clouds to the drag | [`figA2_fractional_contribution.py`](scripts/figA2_fractional_contribution.py) | analytic |
| A.3 | Power-law vs delta cloud models | [`figA3_powerlaw_vs_delta.py`](scripts/figA3_powerlaw_vs_delta.py) | analytic |
| B.1 | 1D stationary outflow tests | [`figB1_1d_tests.py`](scripts/figB1_1d_tests.py) | `data/1d_tests` |
| D.1 | Resolution test | [`figD1_resolution.py`](scripts/figD1_resolution.py) | `data/profiles` |
| E.1 | Hot-phase heating terms | [`fig06_09_E1_E2_mixing.py`](scripts/fig06_09_E1_E2_mixing.py) | `data/profiles` |
| E.2 | Hot-phase temperature and total heating | ″ | `data/profiles` |
| F.1 | Varying *f*<sub>cold</sub>; delta cloud model | [`figF1_fcold_rcl.py`](scripts/figF1_fcold_rcl.py) | `data/profiles` |

The output of each script is written to [`figures/`](figures).

## Repository layout

```
scripts/
  common.py            constants, code units, plot style, snapshot reader,
                       slice sampling and radial binning shared by all scripts
  cloud_model.py       analytic coefficients of the power-law cloud model (Eqs. 18-27)
  make_profiles.py     snapshots -> time-averaged radial profiles (data/profiles)
  fig*.py              one script per figure (or group of related figures)
data/
  profiles/<run>/      time-averaged radial profiles of each run (.npy)
  schneider2020/       curves digitised from Schneider et al. (2020), Figs. 6-8
  1d_tests/            1D test runs and semi-analytic solutions (Appendix B)
  fid_mcrit_histogram.npz   m_crit histogram of the fid run (Fig. A.1, top)
figures/               figures produced by the scripts
```

## Usage

Requires Python ≥ 3.9 and

```bash
pip install -r requirements.txt
```

The figures built from the processed data or from analytic expressions can be
made directly:

```bash
python scripts/fig06_09_E1_E2_mixing.py
```

The same holds for `figA1_rinv_coefficients.py`, `figA2_fractional_contribution.py`,
`figA3_powerlaw_vs_delta.py`, `figB1_1d_tests.py`, `figD1_resolution.py` and
`figF1_fcold_rcl.py`.

Figs. 2–5 (and the m<sub>crit</sub> histogram of Fig. A.1) need the AREPO snapshots.
The scripts expect one folder per run, containing `output/snap_XXX.hdf5`, inside a
root directory given by `--sim-root` or the `M82_SIM_DIR` environment variable:

```bash
export M82_SIM_DIR=/path/to/Simulations
python scripts/fig02_density_slices.py
python scripts/fig03_04_05_schneider_comparison.py
python scripts/make_profiles.py fmix0.1     # regenerate data/profiles/fmix0.1
```

The snapshots are not included in this repository because of their size.

## Analysis conventions

* All time-averaged quantities use snapshots 119–128, i.e. the interval
  *T* = (27, 30) Myr after the onset of SN injection.
* Velocity and density profiles are medians (shaded: 16th–84th percentiles)
  over 0.3 kpc radial bins, within a cone of semi-opening angle 30° around the
  outflow axis, sampled on an edge-on slice through the galaxy centre.
* Mass outflow rates are Ṁ(r) = Σᵢ mᵢ vᵢ / Δr over all cells in the radial bin
  and cone, averaged over the snapshots (shaded: standard deviation in time).
* The cold-phase density in Fig. 4 is the density of the individual sub-grid
  clouds, ρ<sub>cl</sub> = χ ρ<sub>hot</sub> (Eq. 11).

### Run names

Runs are named as in Table 2 of the paper; [`common.RUN_DIRS`](scripts/common.py)
maps them to the simulation folder names:

| Paper | Folder | | Paper | Folder |
|---|---|---|---|---|
| `fid` | `fid` | | `Vc0.001` | `high_res_rainer` |
| `NoSources` | `15J_NoSources` | | `Vc0.01` | `HR1` |
| `NoDrag_fmix0.03` | `fc0.03_nd` | | `Vc0.03` | `HR_Vc0.03` |
| `fmix0.03` | `fc0.03` | | `Vc0.3` | `HR_Vc0.3` |
| `fmix0.1` | `fc0.1` | | `Vc1` | `HR_Vc1` |
| `fmix0.3` | `fc0.3` | | `fcold0.1` | `fcld_l` |
| `fmix1` | `fc1` | | `fcold0.9` | `fcld_h` |
| `rcl0.01` | `rcl0.01` | | `rcl0.6` | `rcl0.6` |

## Citation

If you use this code or data, please cite:

```bibtex
@ARTICLE{Bollati2026,
       author = {{Bollati}, F. and {Weinberger}, R.},
        title = "{Modeling multiphase galactic outflows: A multifluid moving-mesh approach}",
      journal = {Astronomy \& Astrophysics},
         year = 2026,
       volume = {710},
        pages = {A88},
          doi = {10.1051/0004-6361/202556959}
}
```

The curves in `data/schneider2020/` are digitised from Schneider, E. E.,
Ostriker, E. C., Robertson, B. E., & Thompson, T. A. 2020, ApJ, 895, 43;
please cite that work as well if you use them.
