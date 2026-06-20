"""Plotting utilities for Quarles et al. (2018) stability-map files."""

from __future__ import annotations

from pathlib import Path
import re
import tarfile
from typing import Iterable

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from scipy.interpolate import griddata


def read_maxecc_file(path_or_file) -> np.ndarray:
    """Read a ``MaxEcc`` text file from a path or file-like object."""
    data = np.genfromtxt(path_or_file, delimiter=',', comments='#')
    if data.ndim != 2 or data.shape[1] < 5:
        raise ValueError('Expected columns: a_0, M_0/e_0, emax, emin, tcol')
    return data


def find_tar_member(tar_path: str | Path, mu: float, e_bin: float) -> str:
    """Return the member name for a ``MaxEcc_[mu, e_bin].txt`` file in the Zenodo tarball."""
    fname = f'MaxEcc_[{mu:0.3f},{e_bin:0.3f}].txt'
    with tarfile.open(tar_path, 'r:gz') as tar:
        names = tar.getnames()
    if fname in names:
        return fname
    # Some older sample files used one decimal for labels such as 0.0.
    candidates = [name for name in names if name.endswith(fname)]
    if candidates:
        return candidates[0]
    raise FileNotFoundError(f'Could not find {fname} in {tar_path}')


def read_maxecc_from_tar(tar_path: str | Path, mu: float, e_bin: float) -> np.ndarray:
    """Read one ``MaxEcc`` member from ``MaxEcc.tar.gz`` without extracting it."""
    member = find_tar_member(tar_path, mu, e_bin)
    with tarfile.open(tar_path, 'r:gz') as tar:
        with tar.extractfile(member) as handle:
            return read_maxecc_file(handle)


def critical_ac_from_phase_grid(data: np.ndarray, tscale: float = 1e5, tolerance: float = 1e-8) -> float:
    """Return the smallest semimajor axis stable for every sampled phase.

    Parameters
    ----------
    data
        Columns are ``a_0``, phase, ``emax``, ``emin``, and ``tcol``.
    tscale
        Full integration time. Rows with ``tcol == tscale`` are treated as stable.
    """
    a_values = np.array(sorted(np.unique(np.round(data[:, 0], 6))))
    candidate_values = []
    for a0 in a_values:
        rows = np.where(np.isclose(data[:, 0], a0, atol=tolerance))[0]
        if len(rows) == 0:
            continue
        if np.all(np.isclose(data[rows, -1], tscale, atol=tolerance)):
            candidate_values.append(a0)
    if not candidate_values:
        return np.nan
    return float(np.min(candidate_values))


def _grid_for_phase_plot(data: np.ndarray, value: np.ndarray, x_step: float | None = None, y_step: float = 0.01, method: str = 'nearest'):
    x = data[:, 1]
    y = data[:, 0]
    if x_step is None:
        ux = np.unique(x)
        x_step = np.median(np.diff(np.sort(ux))) if len(ux) > 1 else 1.0
    xi = np.arange(np.nanmin(x), np.nanmax(x) + 0.5*x_step, x_step)
    yi = np.arange(np.nanmin(y), np.nanmax(y) + 0.5*y_step, y_step)
    zi = griddata((x, y), value, (xi[None, :], yi[:, None]), method=method)
    return xi, yi, zi


def plot_phase_map(data: np.ndarray, output: str | Path | None = None, tscale: float = 1e5, title: str | None = None, mode: str = 'emax', dpi: int = 300):
    """Plot a single semimajor-axis/phase map.

    ``mode='emax'`` plots ``log10(e_max)`` with unstable cells marked white.
    ``mode='time'`` plots ``log10(|t_col|)``.
    ``mode='delta_e'`` plots ``e_max-e_min``.
    """
    stable = np.isclose(data[:, -1], tscale)
    ac = critical_ac_from_phase_grid(data, tscale=tscale)

    if mode == 'emax':
        values = np.array(data[:, 2], dtype=float)
        values[~stable] = np.nan
        values = np.log10(np.clip(values, 1e-4, None))
        label = r'$\log_{10}(e_{\max})$'
        vmin, vmax = -3.0, 0.0
        cmap = cm.gnuplot_r.copy()
        cmap.set_bad('white')
    elif mode == 'time':
        values = np.log10(np.clip(np.abs(data[:, -1]), 1e-8, None))
        label = r'$\log_{10}|t_{\rm col}|$'
        vmin, vmax = 0.0, np.nanmax(values)
        cmap = cm.gist_rainbow.copy()
    elif mode == 'delta_e':
        values = np.clip(data[:, 2] - data[:, 3], 1e-4, None)
        values[~stable] = np.nan
        values = np.log10(values)
        label = r'$\log_{10}(\Delta e)$'
        vmin, vmax = -3.0, 0.0
        cmap = cm.gnuplot_r.copy()
        cmap.set_bad('white')
    else:
        raise ValueError("mode must be 'emax', 'time', or 'delta_e'")

    xi, yi, zi = _grid_for_phase_plot(data, values)
    fig, ax = plt.subplots(figsize=(8, 7), dpi=dpi)
    mesh = ax.pcolormesh(xi, yi, zi, shading='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    if np.isfinite(ac):
        ax.axhline(ac, color='k', lw=2)
        ax.text(0.02, 0.03, fr'$a_c={ac:.3f}$', transform=ax.transAxes, weight='bold')
    if title:
        ax.set_title(title)
    ax.set_xlabel('Mean anomaly (deg)')
    ax.set_ylabel(r'$a_p/a_{\rm bin}$')
    ax.set_ylim(np.nanmin(data[:, 0]), np.nanmax(data[:, 0]))
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(label)
    fig.tight_layout()
    if output is not None:
        fig.savefig(output, dpi=dpi, bbox_inches='tight')
    return fig, ax


def plot_archive_two_panel(data: np.ndarray, output: str | Path | None = None, dpi: int = 300):
    """Recreate the later two-panel Zenodo diagnostic: lifetime and eccentricity range."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 14), dpi=dpi, sharex=True)
    x = data[:, 1]
    y = data[:, 0]
    xi = np.arange(np.nanmin(x), np.nanmax(x) + 1, 1.0)
    yi = np.arange(np.nanmin(y), np.nanmax(y) + 0.01, 0.01)

    zt = griddata((x, y), np.log10(np.clip(np.abs(data[:, -1]), 1e-8, None)), (xi[None, :], yi[:, None]), method='linear')
    im1 = ax1.pcolormesh(xi, yi, zt, shading='auto', cmap=cm.gist_rainbow)
    ax1.set_ylabel(r'$a_p/a_{\rm bin}$')
    cbar1 = fig.colorbar(im1, ax=ax1)
    cbar1.set_label(r'$\log_{10}|t_{\rm col}|$')

    ze = griddata((x, y), data[:, 2] - data[:, 3], (xi[None, :], yi[:, None]), method='linear')
    im2 = ax2.pcolormesh(xi, yi, ze, shading='auto', cmap=cm.gist_rainbow, vmin=0, vmax=1)
    ax2.set_xlabel('Mean anomaly (deg)')
    ax2.set_ylabel(r'$a_p/a_{\rm bin}$')
    cbar2 = fig.colorbar(im2, ax=ax2)
    cbar2.set_label(r'$\Delta e$')
    fig.tight_layout()
    if output is not None:
        fig.savefig(output, dpi=dpi, bbox_inches='tight')
    return fig, (ax1, ax2)
