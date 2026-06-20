"""Paper-style figure helpers for Quarles et al. (2018).

These routines are intentionally close to the published figure layouts. They
are meant for the example notebooks, not as a replacement for the original
analysis scripts.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from scipy.interpolate import griddata

from cbp_plotting import read_maxecc_file, read_maxecc_from_tar, critical_ac_from_phase_grid

MJUP_MSUN = 9.54e-4

KEPLER_SYSTEMS = [
    {"name": "16", "label": "Kepler-16", "mA": 0.6897, "mB": 0.20255, "mu": 0.2270, "abin": 0.22431, "ebin": 0.15944, "omega_bin": 263.464, "ma_bin": 18.888, "ap": 0.7048, "ep": 0.0069, "ac": 0.6050, "beta": 2.4610, "mp": 0.333*MJUP_MSUN},
    {"name": "34", "label": "Kepler-34", "mA": 1.0479, "mB": 1.0208, "mu": 0.4934, "abin": 0.22882, "ebin": 0.52087, "omega_bin": 71.437, "ma_bin": 228.760, "ap": 1.0896, "ep": 0.18, "ac": 0.8118, "beta": 7.1703, "mp": 0.22*MJUP_MSUN},
    {"name": "35", "label": "Kepler-35", "mA": 0.8877, "mB": 0.8094, "mu": 0.4769, "abin": 0.17617, "ebin": 0.1421, "omega_bin": 89.1784, "ma_bin": 2.9021, "ap": 0.6035, "ep": 0.042, "ac": 0.4795, "beta": 6.3175, "mp": 0.127*MJUP_MSUN},
    {"name": "38", "label": "Kepler-38", "mA": 0.949, "mB": 0.249, "mu": 0.208, "abin": 0.1469, "ebin": 0.1032, "omega_bin": 268.68, "ma_bin": 181.32, "ap": 0.4644, "ep": 0.032, "ac": 0.4328, "beta": 1.1968, "mp": 0.38*MJUP_MSUN},
    {"name": "47", "label": "Kepler-47", "mA": 0.957, "mB": 0.342, "mu": 0.263, "abin": 0.08145, "ebin": 0.0288, "omega_bin": 226.253, "ma_bin": 310.818, "ap": 0.2956, "ep": 0.035, "ac": 0.1848, "beta": 13.519, "mp": 0.1*MJUP_MSUN},
    {"name": "64", "label": "Kepler-64", "mA": 1.528, "mB": 0.408, "mu": 0.211, "abin": 0.1744, "ebin": 0.2117, "omega_bin": 219.7504, "ma_bin": 251.58, "ap": 0.6340, "ep": 0.054, "ac": 0.5368, "beta": 2.9697, "mp": 0.211*MJUP_MSUN},
    {"name": "413", "label": "Kepler-413", "mA": 0.820, "mB": 0.5423, "mu": 0.398, "abin": 0.10148, "ebin": 0.0365, "omega_bin": 279.54, "ma_bin": 169.5328, "ap": 0.3530, "ep": 0.12, "ac": 0.2389, "beta": 8.3487, "mp": 0.211*MJUP_MSUN},
    {"name": "453", "label": "Kepler-453", "mA": 0.944, "mB": 0.1951, "mu": 0.171, "abin": 0.185319, "ebin": 0.0524, "omega_bin": 263.05, "ma_bin": 187.705, "ap": 0.7903, "ep": 0.038, "ac": 0.4184, "beta": 20.152, "mp": 0.03*MJUP_MSUN},
    {"name": "1647", "label": "Kepler-1647", "mA": 1.2207, "mB": 0.9678, "mu": 0.4422, "abin": 0.1276, "ebin": 0.1602, "omega_bin": 300.5442, "ma_bin": 139.0749, "ap": 2.7200, "ep": np.nan, "ac": 0.3497, "beta": 20.275, "mp": 1.52*MJUP_MSUN},
]

SYSTEM_BY_NAME = {row["name"]: row for row in KEPLER_SYSTEMS}

FIG8_ORDER = [["16", "47"], ["34", "64"], ["35", "413"], ["38", "453"]]
FIG8_NAMES = ["16", "34", "35", "38", "47", "64", "413", "453"]
FIG8_OBSERVED = {
    "16": (0.70, 0.0069),
    "34": (1.09, 0.18),
    "35": (0.60, 0.042),
    "38": (0.46, 0.032),
    "47": (0.27, 0.035),
    "64": (0.63, 0.054),
    "413": (0.36, 0.12),
    "453": (0.79, 0.038),
}
FIG8_MU = {
    "16": 0.20255/(0.20255 + 0.68970),
    "34": 1.0208/(1.0208 + 1.0479),
    "35": 0.8094/(0.8094 + 0.8877),
    "38": 0.249/(0.249 + 0.949),
    "47": 0.342/(0.342 + 0.957),
    "64": 0.408/(0.408 + 1.528),
    "413": 0.5423/(0.5423 + 0.82),
    "453": 0.1951/(0.1951 + 0.944),
}
FIG8_EBIN = {
    "16": 0.15944,
    "34": 0.52087,
    "35": 0.1421,
    "38": 0.1032,
    "47": 0.0288,
    "64": 0.2117,
    "413": 0.0365,
    "453": 0.0524,
}
FIG8_ABIN = {
    "16": 0.22431,
    "34": 0.22882,
    "35": 0.17617,
    "38": 0.1649,
    "47": 0.08145,
    "64": 0.1744,
    "413": 0.10148,
    "453": 0.185319,
}
FIG8_INTERP_AC_AU = {
    "16": 0.60,
    "34": 0.81,
    "35": 0.48,
    "38": 0.43,
    "47": 0.18,
    "64": 0.52,
    "413": 0.24,
    "453": 0.42,
}

FIT_HW99 = [1.60, 5.10, -2.22, 4.12, -4.27, -5.09, 4.61]
FIT_Q18_1 = [1.48, 3.92, -1.41, 5.14, 0.33, -7.95, -4.89]


def stability_fit(mu: float | np.ndarray, ebin: float | np.ndarray, coeffs: Iterable[float]) -> np.ndarray:
    c1, c2, c3, c4, c5, c6, c7 = coeffs
    return c1 + c2*ebin + c3*ebin**2 + c4*mu + c5*mu*ebin + c6*mu**2 + c7*(mu*ebin)**2


def nearest_archive_grid(system: dict) -> tuple[float, float]:
    mu_grid = round(float(system["mu"]), 2)
    e_grid = round(float(system["ebin"]), 2)
    if mu_grid == 0.0:
        mu_grid = 0.001
    return mu_grid, e_grid


def folded_angle(angle_deg: float) -> float:
    angle = float(angle_deg) % 360.0
    return min(angle, 360.0 - angle)


def plot_figure5(a_crit_path: str | Path = "a_crit.txt", output: str | Path | None = None, dpi: int = 300):
    """Recreate paper Figure 5 from ``a_crit.txt``."""
    data = np.genfromtxt(a_crit_path, delimiter=',', comments='#')
    x = data[:, 0].copy()
    x[np.isclose(x, 0.001)] = 0.0
    y = data[:, 1]
    z = data[:, 2]
    xi = np.linspace(0.0, 0.5, 51)
    yi = np.linspace(0.0, 0.8, 81)
    zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='linear', fill_value=np.nan)

    cmap = cm.nipy_spectral.copy()
    cmap.set_bad('white')
    fig, ax = plt.subplots(figsize=(12, 6.75), dpi=dpi)
    mesh = ax.pcolormesh(xi, yi, zi, shading='auto', cmap=cmap, vmin=1.3, vmax=4.5)
    for system in KEPLER_SYSTEMS:
        ax.plot(system['mu'], system['ebin'], '.', color='white', ms=10)
        ax.text(system['mu'], system['ebin'] + 0.02, system['name'], color='white', fontsize='medium', ha='center')
    ax.set_xlabel(r'$\mu$')
    ax.set_ylabel(r'$e_{\rm bin}$')
    ax.set_xlim(0.0, 0.5)
    ax.set_ylim(0.0, 0.8)
    ax.tick_params(axis='both', direction='out', length=4.0, width=1.5)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r'$a_c$')
    fig.tight_layout()
    if output is not None:
        fig.savefig(output, dpi=dpi, bbox_inches='tight')
    return fig, ax


def plot_figure6(output: str | Path | None = None, dpi: int = 300):
    """Recreate the Figure 6 dynamical-spacing schematic."""
    fig, ax = plt.subplots(figsize=(7.6, 7.6), dpi=dpi)
    theta = np.linspace(0.0, 0.5*np.pi, 500)
    label_angles = {"38": 40, "16": 50, "64": 25, "35": 45, "413": 55, "47": 60, "453": 35, "1647": 60}

    for system in KEPLER_SYSTEMS:
        beta = system['beta']
        color = 'red' if beta < 7.0 else 'black'
        linestyle = '-'
        if system['name'] == '453':
            linestyle = '--'
        if system['name'] == '1647':
            linestyle = '-'
        ax.plot(beta*np.cos(theta), beta*np.sin(theta), color=color, linestyle=linestyle, lw=1.5)
        angle = np.deg2rad(label_angles.get(system['name'], 45.0))
        ax.text(beta*np.cos(angle), beta*np.sin(angle), system['name'], color=color, fontsize=10, ha='center', va='center', rotation=np.rad2deg(angle)-35)

    beta = 7.0
    ax.plot(beta*np.cos(theta), beta*np.sin(theta), color='red', lw=1.5)
    ax.text(beta*np.cos(np.deg2rad(78)), beta*np.sin(np.deg2rad(78)), '7', color='red', fontsize=10, ha='center', va='center')
    ax.text(0.98, 0.97, r'$\beta_c = \frac{a_p-a_c}{R_{H,m}}$', transform=ax.transAxes, ha='right', va='top')
    ax.set_xlim(0.0, 20.3)
    ax.set_ylim(0.0, 20.3)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(r'$\beta_c$ (Mutual Hill Radii)')
    ax.set_ylabel(r'$\beta_c$ (Mutual Hill Radii)')
    ax.set_xticks(np.arange(0, 21, 5))
    ax.set_yticks(np.arange(0, 21, 5))
    ax.tick_params(axis='both', which='both', direction='out', top=True, right=True, length=4.0, width=1.2)
    fig.tight_layout()
    if output is not None:
        fig.savefig(output, dpi=dpi, bbox_inches='tight')
    return fig, ax


def _map_values_for_emax(data: np.ndarray, tscale: float = 1e5):
    stable = np.isclose(data[:, -1], tscale)
    values = np.asarray(data[:, 2], dtype=float).copy()
    values[~stable] = np.nan
    return np.log10(np.clip(values, 1e-3, None))


def _grid_xy(data: np.ndarray, values: np.ndarray, x_step: float, y_step: float, method: str = 'nearest'):
    x = data[:, 1]
    y = data[:, 0]
    xi = np.arange(np.nanmin(x), np.nanmax(x) + 0.5*x_step, x_step)
    yi = np.arange(np.nanmin(y), np.nanmax(y) + 0.5*y_step, y_step)
    zi = griddata((x, y), values, (xi[None, :], yi[:, None]), method=method)
    return xi, yi, zi


def plot_figure7(tar_path: str | Path = "MaxEcc.tar.gz", output: str | Path | None = None, dpi: int = 300):
    """Recreate paper Figure 7 from the Zenodo ``MaxEcc.tar.gz`` archive."""
    tar_path = Path(tar_path)
    if not tar_path.exists():
        raise FileNotFoundError(f"{tar_path} was not found. Download MaxEcc.tar.gz from the Zenodo archive and place it in the repository root.")

    order = [["16", "38", "413"], ["34", "47", "453"], ["35", "64", "1647"]]
    cmap = cm.gnuplot_r.copy()
    cmap.set_bad('white')
    norm = colors.Normalize(-3.0, 0.0)
    fig, axes = plt.subplots(3, 3, figsize=(12, 8), dpi=dpi, sharex=True, sharey=True)
    last_mesh = None

    for i, row in enumerate(order):
        for j, name in enumerate(row):
            ax = axes[i, j]
            system = SYSTEM_BY_NAME[name]
            mu_grid, e_grid = nearest_archive_grid(system)
            data = read_maxecc_from_tar(tar_path, mu_grid, e_grid)
            values = _map_values_for_emax(data, tscale=1e5)
            xi, yi, zi = _grid_xy(data, values, x_step=2.0, y_step=0.01, method='nearest')
            last_mesh = ax.pcolormesh(xi, yi, zi, shading='auto', cmap=cmap, norm=norm)
            ac = critical_ac_from_phase_grid(data, tscale=1e5)
            if np.isfinite(ac):
                ax.axhline(ac, color='c', lw=2)
                ax.text(0.02, 0.04, fr'$a_c={ac:.2f}$', color='c', transform=ax.transAxes, fontsize=9, weight='bold')
            ydot = system['ap']/system['abin']
            xdot = folded_angle(system['ma_bin'])
            if ydot <= 5.0:
                ax.plot(xdot, ydot, '.', color='limegreen', ms=8)
            ax.text(0.02, 0.90, f"({mu_grid:.2f},{e_grid:.2f})", transform=ax.transAxes, color='black', fontsize=10, weight='bold')
            ax.text(0.02, 0.78, system['label'], transform=ax.transAxes, color='limegreen', fontsize=10, weight='bold')
            ax.set_xlim(0, 180)
            ax.set_ylim(1, 5)
            ax.set_xticks(np.arange(0, 181, 30))
            ax.set_yticks(np.arange(1, 6, 1))
            ax.tick_params(axis='both', direction='out', length=3.0, width=1.0)

    fig.text(0.5, 0.04, 'Mean Anomaly (deg.)', ha='center')
    fig.text(0.04, 0.5, r'$a_p/a_{\rm bin}$', va='center', rotation='vertical')
    axes[0, 2].text(1.0, 1.05, r'$\lambda_{\rm bin}=0^\circ$', transform=axes[0, 2].transAxes, ha='right', fontsize=10, weight='bold')
    fig.subplots_adjust(left=0.09, right=0.86, bottom=0.10, top=0.95, hspace=0.08, wspace=0.08)
    cax = fig.add_axes([0.88, 0.15, 0.018, 0.75])
    cbar = fig.colorbar(last_mesh, cax=cax)
    cbar.set_label(r'$\log_{10}(e_{\max})$')
    if output is not None:
        fig.savefig(output, dpi=dpi, bbox_inches='tight')
    return fig, axes


def _grid_ap_ep(data: np.ndarray, values: np.ndarray, x_step: float = 0.001, y_step: float = 0.01):
    x = data[:, 0]
    y = data[:, 1]
    xi = np.arange(np.nanmin(x), np.nanmax(x) + 0.5*x_step, x_step)
    yi = np.arange(np.nanmin(y), np.nanmax(y) + 0.5*y_step, y_step)
    zi = griddata((x, y), values, (xi[None, :], yi[:, None]), method='nearest')
    return xi, yi, zi


def plot_figure8(data_dir: str | Path = "plot_figures/Fig8_data", output: str | Path | None = None, dpi: int = 300, two_planet_tmin: float = 1e7):
    """Recreate paper Figure 8 from the included Figure 8 data files.

    This function intentionally follows the original ``plot_Fig8_single.py``
    conventions used for the published figure. In particular, the stability
    curves use the original Figure 8 constants, the 0.8 pericenter factor, and
    the rounded interpolation values plotted in the paper figure rather than
    the more precise Table 4 values. The gray two-planet markers are restricted
    to interior-planet semimajor axes inside the observed planet pericenter,
    ``q_p = a_p(1-e_p)``, matching the setup described for the packing tests.
    """
    data_dir = Path(data_dir)
    if not data_dir.exists():
        raise FileNotFoundError(f"{data_dir} was not found.")

    cmap = cm.gnuplot_r.copy()
    cmap.set_under('gray')
    cmap.set_over('white')
    cmap.set_bad('white')
    norm = colors.Normalize(-3.0, 0.0)
    fig, axes = plt.subplots(4, 2, figsize=(11.0, 7.8), dpi=dpi, sharex=True, sharey=True)
    legend_handles = []
    last_mesh = None
    fact = 1.0 - 0.2

    for i, row in enumerate(FIG8_ORDER):
        for j, name in enumerate(row):
            ax = axes[i, j]
            single_path = data_dir / f"MaxEcc_{name}.txt"
            data = read_maxecc_file(single_path)
            data = np.asarray(data, dtype=float)
            unstab = data[:, -1] < 1e5
            data_plot = data.copy()
            data_plot[unstab, 2] = 3.0
            data_plot[unstab, 3] = 3.0
            x = data_plot[:, 0]
            y = data_plot[:, 1]
            z = data_plot[:, 2] - data_plot[:, 3]
            z[unstab] = 3.0
            xmin, xmax = np.nanmin(x), np.nanmax(x)
            nx = int((xmax - xmin)/0.001) + 1
            xi = np.linspace(xmin, xmax, nx)
            yi = np.linspace(0.0, 0.5, 51)
            zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='nearest')
            last_mesh = ax.pcolormesh(xi, yi, np.log10(np.clip(zi, 1e-12, None)), cmap=cmap, norm=norm, shading='auto', zorder=2)

            two_path = data_dir / f"MaxEcc_2pl_{name}.txt"
            if two_path.exists():
                two = read_maxecc_file(two_path)
                two = np.asarray(two, dtype=float)
                obs_a, obs_e = FIG8_OBSERVED[name]
                q_obs = obs_a*(1.0 - obs_e)
                two_stable = (two[:, -1] >= two_planet_tmin) & (two[:, 0] <= q_obs + 1e-12)
                if np.any(two_stable):
                    ax.scatter(two[two_stable, 0], two[two_stable, 1], marker='s', s=8, color='0.65', edgecolors='none', zorder=4)

            mu = FIG8_MU[name]
            ebin = FIG8_EBIN[name]
            abin = FIG8_ABIN[name]
            ac_fit1 = stability_fit(mu, ebin, FIT_Q18_1)*abin
            ac_hw99 = stability_fit(mu, ebin, FIT_HW99)*abin
            ac_interp = FIG8_INTERP_AC_AU[name]
            y_fit1 = fact*(1.0 - ac_fit1/xi)
            y_interp = fact*(1.0 - ac_interp/xi)
            y_hw99 = fact*(1.0 - ac_hw99/xi)
            h1, = ax.plot(xi, y_fit1, color='c', lw=2.0, label='Fit 1', zorder=5)
            h2, = ax.plot(xi, y_interp, color='y', lw=2.0, linestyle='-.', label='Interpolation', zorder=5)
            h3, = ax.plot(xi, y_hw99, color='violet', lw=2.0, linestyle='--', label='HW99', zorder=5)
            if not legend_handles:
                legend_handles = [h1, h2, h3]

            obs_a, obs_e = FIG8_OBSERVED[name]
            ax.plot(obs_a, obs_e, '.', color='limegreen', ms=10, zorder=6)
            ax.text(0.96, 0.07, f"Kepler-{name}", transform=ax.transAxes, color='limegreen', fontsize=10, weight='bold', ha='right')
            ax.set_xlim(0.15, 1.5)
            ax.set_ylim(0.0, 0.5)
            ax.set_xticks([0.2 + k*0.1 for k in range(14)])
            ax.set_yticks([k*0.1 for k in range(6)])
            ax.tick_params(axis='both', direction='out', length=3.0, width=1.0)

    fig.legend(handles=legend_handles, loc='upper center', ncol=3, frameon=False, bbox_to_anchor=(0.46, 0.99), handletextpad=0.25)
    fig.text(0.5, 0.04, r'$a_p$ (AU)', ha='center')
    fig.text(0.04, 0.5, r'$e_p$', va='center', rotation='vertical')
    fig.subplots_adjust(left=0.09, right=0.86, bottom=0.09, top=0.94, hspace=0.03, wspace=0.03)
    cax = fig.add_axes([0.88, 0.13, 0.018, 0.74])
    cbar = fig.colorbar(last_mesh, cax=cax)
    cbar.set_label(r'$\log_{10}[\Delta e]$')
    if output is not None:
        fig.savefig(output, dpi=dpi, bbox_inches='tight')
    return fig, axes

