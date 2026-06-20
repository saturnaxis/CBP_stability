"""Utilities for the Quarles et al. (2018) circumbinary stability grid.

The published data file, ``a_crit.txt``, tabulates the critical semimajor-axis
ratio ``a_c = a_p/a_bin`` as a function of binary mass ratio ``mu`` and binary
eccentricity ``e_bin``.  These helpers provide a Python 3 interface to that grid
without relying on the deprecated ``scipy.interpolate.interp2d`` function.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.interpolate import RegularGridInterpolator


@dataclass(frozen=True)
class StabilityGrid:
    """Container for the regular ``a_c(mu, e_bin)`` interpolation grid."""

    mu_values: np.ndarray
    eb_values: np.ndarray
    ac_values: np.ndarray
    interpolator: RegularGridInterpolator


def q_to_mu(q: float) -> float:
    """Convert observational mass quotient ``q = M_B/M_A`` to ``mu``."""
    return float(q) / (1.0 + float(q))


def mu_from_masses(m1: float, m2: float) -> float:
    """Return the smaller-star mass fraction used by Quarles et al. (2018)."""
    m1 = float(m1)
    m2 = float(m2)
    return min(m1, m2) / (m1 + m2)


def _repo_root() -> Path:
    return Path(__file__).resolve().parent


def _normalize_mu(mu: float) -> float:
    """Map the formal low-mass limit to the tabulated ``mu=0.001`` row."""
    mu = float(mu)
    if mu < 0.001:
        return 0.001
    return mu


def load_acrit(path: str | Path | None = None) -> np.ndarray:
    """Load the ``a_crit.txt`` table.

    Returns
    -------
    numpy.ndarray
        Columns are ``mu``, ``e_bin``, and ``a_c``.
    """
    if path is None:
        path = _repo_root() / 'a_crit.txt'
    data = np.genfromtxt(path, delimiter=',', comments='#')
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f'Could not read a three-column a_crit table from {path!s}')
    return data[:, :3]


def build_stability_grid(path: str | Path | None = None) -> StabilityGrid:
    """Build a regular-grid interpolator from ``a_crit.txt``.

    The original table is unordered and contains one duplicated grid point.  Duplicate
    entries are averaged before the regular grid is constructed.
    """
    data = load_acrit(path)
    mu_values = np.unique(data[:, 0])
    eb_values = np.unique(data[:, 1])
    ac_values = np.full((len(mu_values), len(eb_values)), np.nan, dtype=float)

    for i, mu in enumerate(mu_values):
        for j, eb in enumerate(eb_values):
            mask = np.isclose(data[:, 0], mu) & np.isclose(data[:, 1], eb)
            if np.any(mask):
                ac_values[i, j] = np.nanmean(data[mask, 2])

    if np.any(np.isnan(ac_values)):
        missing = np.argwhere(np.isnan(ac_values))
        raise ValueError(f'a_crit grid has missing values at indices {missing[:5].tolist()}')

    interpolator = RegularGridInterpolator((mu_values, eb_values), ac_values, bounds_error=False, fill_value=np.nan)
    return StabilityGrid(mu_values=mu_values, eb_values=eb_values, ac_values=ac_values, interpolator=interpolator)


def get_ac(mu: float, e_bin: float, grid: StabilityGrid | None = None) -> float:
    """Interpolate the critical semimajor-axis ratio ``a_c = a_p/a_bin``."""
    mu = _normalize_mu(mu)
    e_bin = float(e_bin)
    if grid is None:
        grid = build_stability_grid()
    if mu > float(grid.mu_values.max()) or e_bin < float(grid.eb_values.min()) or e_bin > float(grid.eb_values.max()):
        raise ValueError('Requested parameters are outside the tabulated domain: 0 <= mu <= 0.5 and 0 <= e_bin <= 0.8.')
    value = grid.interpolator([[mu, e_bin]])[0]
    if not np.isfinite(value):
        raise ValueError('Interpolation returned NaN. Check that the input is inside the tabulated grid.')
    return float(value)


def critical_period_days(mu: float, e_bin: float, p_bin_days: float, grid: StabilityGrid | None = None) -> float:
    """Return the critical planetary period in days using Kepler's third law."""
    ac = get_ac(mu, e_bin, grid=grid)
    return float(ac**1.5 * float(p_bin_days))


def hw99_ac(mu: float, e_bin: float) -> float:
    """Holman & Wiegert (1999) P-type critical semimajor-axis fit."""
    mu = float(mu)
    e = float(e_bin)
    return 1.60 + 5.10*e - 2.22*e**2 + 4.12*mu - 4.27*mu*e - 5.09*mu**2 + 4.61*mu**2*e**2


def q2018_fit1_ac(mu: float, e_bin: float) -> float:
    """Quarles et al. (2018) polynomial fit using the same terms as HW99."""
    mu = float(mu)
    e = float(e_bin)
    return 1.48 + 3.92*e - 1.41*e**2 + 5.14*mu + 0.33*mu*e - 7.95*mu**2 - 4.89*(mu*e)**2


def q2018_fit2_ac(mu: float, e_bin: float) -> float:
    """Quarles et al. (2018) Hill-motivated polynomial fit."""
    mu = float(mu)
    e = float(e_bin)
    x = mu**(1.0/3.0)
    return 0.93 + 2.67*e - 0.25*e**2 + 3.72*x + 2.25*x*e - 2.72*x**2 - 4.17*x**2*e**2
