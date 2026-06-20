# REBOUND slice examples for Quarles et al. (2018) CBP stability

These scripts are example, slice-wise recreations of the Quarles et al. (2018) workflow using REBOUND. They preserve the original script style where useful: one script per experiment, shuffled grids, job slicing, `MaxEcc_...txt` outputs, `Beta_...txt` outputs, and `a_crit` post-processing.

They are **not** a replacement for the published Mercury6 production archive or the repository `a_crit.txt` table. They are intended for reproducibility checks, teaching, and exploratory extensions.

## What changed after reviewing the additional original scripts

The additional project scripts added three missing pieces to this recreation package:

1. `run_beta_packing_slice.py` and `plot_beta_packing_slice.py` reproduce the beta-spacing two-planet experiment from `Mercury_CBP_beta.py`.
2. `plot_hw_grid_panel.py` reproduces the 4-by-4 panel plotting style from `plot_muEcc.py` and `plot_muEcc_emax.py`.
3. `plot_acrit_summary.py` and `fit_acrit_surface.py` reproduce the population-level `a_crit` surface, the median/range versus binary eccentricity plot, and the coefficient-fitting workflow from `plot_acr.py`, `plot_acritEcc.py`, and `fit_emcee.py`.

I also added `merge_job_outputs.py` because job-split files must be merged before plotting a full slice or measuring `a_c`.

## Quick smoke test

Use a tiny grid first:

```bash
python run_hw_slice.py --mu 0.300 --eb 0.200 --amin 2.0 --amax 2.2 --da 0.1 --dm 30 --orbits 100 --workers 1
python reduce_hw_acrit.py --orbits 100
python plot_hw_slice.py GenRuns_out/MaxEcc_[0.300,0.200].txt --orbits 100
```

## Paper-like single-planet slice

```bash
python run_hw_slice.py --mu 0.300 --eb 0.200 --amin 1.01 --amax 5.0 --da 0.01 --dm 2 --orbits 100000 --workers 16
python reduce_hw_acrit.py --orbits 100000
python plot_hw_slice.py GenRuns_out/MaxEcc_[0.300,0.200].txt --orbits 100000
```

## Split a slice across jobs

```bash
python launch_examples.py --mode hw --njobs 10
python run_hw_slice.py --mu 0.300 --eb 0.200 --jobid 0 --njobs 10 --workers 8
python run_hw_slice.py --mu 0.300 --eb 0.200 --jobid 1 --njobs 10 --workers 8
```

Then merge the job files:

```bash
python merge_job_outputs.py 'GenRuns_out/MaxEcc_[0.300,0.200]_job*.txt'
```

`reduce_hw_acrit.py` can group job-split files automatically, but plotting a single slice is cleaner after merging.

## 4-by-4 panel plots

```bash
python plot_hw_grid_panel.py --indir GenRuns_out --mode emax --orbits 100000
python plot_hw_grid_panel.py --indir GenRuns_out --mode stable --orbits 100000
```

The default panels use `mu = 0.001, 0.100, 0.300, 0.500` and `e_bin = 0.000, 0.100, 0.300, 0.500`.

## Full or partial acrit survey

```bash
python run_acrit_survey_slice.py --mu-min 0.10 --mu-max 0.30 --dmu 0.10 --eb-min 0.0 --eb-max 0.3 --deb 0.1 --orbits 100000 --workers 8
python reduce_hw_acrit.py --orbits 100000 --outfile a_crit_rebound.txt
python plot_acrit_summary.py --acrit a_crit_rebound.txt --plot-kepler
python fit_acrit_surface.py --acrit a_crit_rebound.txt
```

The fitting script writes `fit_outputs/acrit_fit_coefficients.txt`. Add `--run-emcee` only if `emcee` and `corner` are installed.

## Inclined orbit example

```bash
python run_inclination_slice.py --mu 0.300 --eb 0.200 --inc 45 --ac 2.8 --da 0.01 --dm 2 --orbits 500000 --workers 16
```

## Known Kepler CBP host, one test planet

```bash
python run_kepler_ep_slice.py --system Kepler-16b --da 0.01 --de 0.02 --orbits 100000 --workers 16
```

## Two-planet packing slice in semimajor axis and eccentricity

```bash
python run_two_planet_packing_slice.py --system Kepler-16b --da 0.001 --de 0.02 --orbits 100000 --workers 16
```

The default two-planet phasing places the inner planet at `M=0 deg` and the observed outer CBP at `M=180 deg`, matching the fixed-phase spirit of the original two-planet scripts rather than the full mean-anomaly grid used for the single-planet stability maps.

## Two-planet beta-spacing slice

```bash
python run_beta_packing_slice.py --system Kepler-16b --dbeta 0.05 --de 0.02 --orbits 100000 --workers 16
python plot_beta_packing_slice.py CBP_out/Beta_Kepler16b.txt --orbits 100000
```

This is the closest REBOUND analog of `Mercury_CBP_beta.py`: the outer planet is fixed at `ao`, the inner planet is placed using the mutual-Hill beta transformation, and both planets are assigned the same trial eccentricity.

## Output conventions

All numeric outputs are rounded to three decimals at most.

Single-planet files:

```text
#a_0,M_0,emax,emin,tcol
```

Two-dimensional `a-e` files:

```text
#a_0,e_0,emax,emin,tcol
```

Beta-spacing files:

```text
#beta,e_0,emax,emin,tcol
```

A negative `tcol` marks an escape/ejection-style event.
