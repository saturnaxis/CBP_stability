from __future__ import annotations

"""Shared REBOUND utilities for slice-wise recreations of Quarles et al. (2018).

These examples keep the original Mercury-script workflow style where practical:
parameter grids, job slicing, shuffled task order, CSV-like text outputs, and
columns named like the original MaxEcc files. They do not replace the published
Mercury6 archive or the repository a_crit.txt table.
"""

import math
import os
import threading
from dataclasses import dataclass
from multiprocessing.dummy import Pool
from typing import Iterable, List, Sequence, Tuple

import numpy as np
try:
    import rebound
except ImportError:  # plotting/reduction utilities do not require REBOUND
    rebound = None

G = 4.0 * math.pi**2
M_J = 9.54e-4
M_E = 3.0035e-6
DEFAULT_SEED = 21417


@dataclass(frozen=True)
class RunResult:
    a0: float
    m0: float
    e0: float
    emax: float
    emin: float
    tcol: float
    status: str


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def fmt3(value: float) -> str:
    if isinstance(value, (int, np.integer)):
        return f"{int(value)}"
    if isinstance(value, str):
        return value
    if not np.isfinite(value):
        return "nan"
    return f"{float(value):.3f}"


def write_header(path: str, header: str) -> None:
    ensure_dir(os.path.dirname(path) or ".")
    with open(path, "w", encoding="utf-8") as f:
        f.write(header.rstrip() + "\n")


def append_rows(path: str, rows: Sequence[str], lock: threading.Lock | None = None) -> None:
    if lock is None:
        with open(path, "a", encoding="utf-8") as f:
            f.writelines(rows)
        return
    with lock:
        with open(path, "a", encoding="utf-8") as f:
            f.writelines(rows)


def split_job(parameters: Sequence[tuple], jobid: int, njobs: int) -> List[tuple]:
    if njobs <= 0:
        raise ValueError("njobs must be positive")
    n = len(parameters)
    start = jobid * n // njobs
    end = (jobid + 1) * n // njobs
    return list(parameters[start:end])


def shuffled_parameters(parameters: Sequence[tuple], seed: int = DEFAULT_SEED) -> List[tuple]:
    rng = np.random.default_rng(seed)
    idx = np.arange(len(parameters))
    rng.shuffle(idx)
    return [parameters[i] for i in idx]


def binary_com(sim: rebound.Simulation) -> rebound.Particle:
    p0 = sim.particles[0]
    p1 = sim.particles[1]
    mtot = p0.m + p1.m
    com = rebound.Particle(m=mtot)
    com.x = (p0.m * p0.x + p1.m * p1.x) / mtot
    com.y = (p0.m * p0.y + p1.m * p1.y) / mtot
    com.z = (p0.m * p0.z + p1.m * p1.z) / mtot
    com.vx = (p0.m * p0.vx + p1.m * p1.vx) / mtot
    com.vy = (p0.m * p0.vy + p1.m * p1.vy) / mtot
    com.vz = (p0.m * p0.vz + p1.m * p1.vz) / mtot
    return com


def orbital_period(a: float, mtot: float = 1.0) -> float:
    return math.sqrt(a**3 / mtot)


def make_cbp_sim(mu: float, eb: float, a0: float, m0: float, ep: float = 0.0, inc: float = 0.0, planet_mass: float = M_J, binary_phase_deg: float = 0.0, a_bin: float = 1.0) -> rebound.Simulation:
    if rebound is None:
        raise ImportError("REBOUND is required for simulation scripts. Install it with: pip install rebound")
    if not (0.0 < mu <= 0.5):
        raise ValueError("mu should be mB/(mA+mB) and in the range 0 < mu <= 0.5")
    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    sim.G = G
    mA = 1.0 - mu
    mB = mu
    sim.add(m=mA)
    sim.add(m=mB, a=a_bin, e=eb, inc=0.0, Omega=0.0, omega=0.0, f=np.deg2rad(binary_phase_deg))
    sim.move_to_com()
    com = binary_com(sim)
    sim.add(m=planet_mass, a=a0, e=ep, inc=np.deg2rad(inc), Omega=0.0, omega=0.0, M=np.deg2rad(m0), primary=com)
    sim.move_to_com()
    sim.integrator = "whfast"
    peri = max(a0 * (1.0 - ep), 1e-6)
    dt_tp = math.sqrt(peri**3 / 1.0) / 20.0
    sim.dt = 0.48 * min(orbital_period(a_bin, 1.0) / 20.0, dt_tp)
    return sim


def add_outer_cbp(sim: rebound.Simulation, a_outer: float, e_outer: float, mass_outer: float, m_outer_deg: float = 180.0, inc_outer_deg: float = 0.0) -> None:
    com = binary_com(sim)
    sim.add(m=mass_outer, a=a_outer, e=e_outer, inc=np.deg2rad(inc_outer_deg), Omega=0.0, omega=0.0, M=np.deg2rad(m_outer_deg), primary=com)
    sim.move_to_com()


def evaluate_sim(mu: float, eb: float, a0: float, m0: float, ep: float = 0.0, inc: float = 0.0, planet_mass: float = M_J, tscale: float = 1e5, nout: int = 250, binary_phase_deg: float = 0.0, escape_radius: float = 10.0, a_bin: float = 1.0, outer_planet: dict | None = None) -> RunResult:
    sim = make_cbp_sim(mu, eb, a0, m0, ep=ep, inc=inc, planet_mass=planet_mass, binary_phase_deg=binary_phase_deg, a_bin=a_bin)
    if outer_planet is not None:
        add_outer_cbp(sim, outer_planet["a"], outer_planet.get("e", 0.0), outer_planet.get("m", M_J), outer_planet.get("M", 180.0), outer_planet.get("inc", 0.0))
    emax = 0.0
    emin = 1.0
    tcol = float(tscale)
    status = "complete"
    times = np.linspace(0.0, tscale, max(2, int(nout)))
    crossing_radius = a_bin * (1.0 + eb)
    for t in times[1:]:
        try:
            sim.integrate(t, exact_finish_time=0)
        except Exception:
            return RunResult(a0, m0, ep, round(emax, 3), round(emin, 3), round(float(sim.t), 3), "integrator_error")
        p = sim.particles[2]
        com = binary_com(sim)
        dx = p.x - com.x
        dy = p.y - com.y
        dz = p.z - com.z
        r = math.sqrt(dx * dx + dy * dy + dz * dz)
        if r < crossing_radius:
            tcol = float(sim.t)
            status = "crossing"
            break
        if r > escape_radius:
            tcol = -float(sim.t)
            status = "ejected"
            break
        try:
            orb = p.orbit(primary=com)
            if np.isfinite(orb.e):
                emax = max(emax, float(orb.e))
                emin = min(emin, float(orb.e))
        except Exception:
            tcol = float(sim.t)
            status = "orbit_error"
            break
        if outer_planet is not None and len(sim.particles) > 3:
            p3 = sim.particles[3]
            sep = math.sqrt((p.x - p3.x)**2 + (p.y - p3.y)**2 + (p.z - p3.z)**2)
            rh = 0.5 * (a0 + outer_planet["a"]) * ((planet_mass + outer_planet.get("m", M_J)) / 3.0)**(1.0 / 3.0)
            if sep < 0.5 * rh:
                tcol = float(sim.t)
                status = "planet_encounter"
                break
    if emin == 1.0 and emax == 0.0:
        emin = 0.0
    return RunResult(round(a0, 3), round(m0, 3), round(ep, 3), round(emax, 3), round(emin, 3), round(tcol, 3), status)


def result_to_maxecc_row(result: RunResult, include_e0: bool = False) -> str:
    if include_e0:
        return f"{fmt3(result.a0)},{fmt3(result.e0)},{fmt3(result.emax)},{fmt3(result.emin)},{fmt3(result.tcol)}\n"
    return f"{fmt3(result.a0)},{int(round(result.m0)):03d},{fmt3(result.emax)},{fmt3(result.emin)},{fmt3(result.tcol)}\n"


def run_pool(func, parameters: Sequence[tuple], workers: int = 1) -> list:
    if workers <= 1:
        return [func(p) for p in parameters]
    with Pool(processes=workers) as pool:
        return pool.map(func, parameters)


KEPLER_SYSTEMS = {
    "Kepler-16b": {"mu": 0.227, "eb": 0.159, "m_total": 0.89225, "a_bin": 0.224, "a_planet": 0.705, "m_planet": 105.837 * M_E, "e_planet": 0.007},
    "Kepler-34b": {"mu": 0.493, "eb": 0.521, "m_total": 2.06870, "a_bin": 0.229, "a_planet": 1.090, "m_planet": 69.923 * M_E, "e_planet": 0.182},
    "Kepler-35b": {"mu": 0.477, "eb": 0.142, "m_total": 1.69710, "a_bin": 0.176, "a_planet": 0.604, "m_planet": 40.364 * M_E, "e_planet": 0.042},
    "Kepler-38b": {"mu": 0.208, "eb": 0.103, "m_total": 1.19800, "a_bin": 0.147, "a_planet": 0.464, "m_planet": 122.000 * M_E, "e_planet": 0.032},
    "Kepler-47b": {"mu": 0.263, "eb": 0.029, "m_total": 1.29900, "a_bin": 0.081, "a_planet": 0.296, "m_planet": 25.770 * M_E, "e_planet": 0.021},
    "Kepler-64b": {"mu": 0.211, "eb": 0.212, "m_total": 1.93600, "a_bin": 0.174, "a_planet": 0.634, "m_planet": 69.000 * M_E, "e_planet": 0.054},
    "Kepler-413b": {"mu": 0.398, "eb": 0.037, "m_total": 1.36230, "a_bin": 0.101, "a_planet": 0.353, "m_planet": 67.000 * M_E, "e_planet": 0.118},
    "Kepler-453b": {"mu": 0.171, "eb": 0.052, "m_total": 1.13910, "a_bin": 0.185, "a_planet": 0.790, "m_planet": 16.200 * M_E, "e_planet": 0.035},
    "Kepler-1647b": {"mu": 0.442, "eb": 0.160, "m_total": 2.18850, "a_bin": 0.128, "a_planet": 2.720, "m_planet": 483.000 * M_E, "e_planet": 0.058},
}


def normalize_mass_to_binary(mass_msun: float, sysdata: dict) -> float:
    """Return planet mass in the normalized units where the binary total mass is 1."""
    return float(mass_msun) / float(sysdata.get("m_total", 1.0))
