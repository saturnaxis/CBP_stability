#!/usr/bin/env python3
"""Two-planet packing slice for a Kepler CBP system.

This follows the intent of the original CBP two-planet scripts: keep the observed
outer CBP fixed, insert an equal- or user-selected-mass inner planet, and sweep
inner-planet semimajor axis and eccentricity. The default phasing puts the inner
planet at M=0 deg and the observed outer planet at M=180 deg.
"""

import argparse
import time

import numpy as np

from cbp_common import KEPLER_SYSTEMS, M_E, append_rows, ensure_dir, evaluate_sim, normalize_mass_to_binary, result_to_maxecc_row, run_pool, shuffled_parameters, split_job, write_header


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", type=str, default="Kepler-16b", choices=sorted(KEPLER_SYSTEMS))
    parser.add_argument("--amin", type=float, default=None, help="Minimum inner-planet semimajor axis in AU.")
    parser.add_argument("--amax", type=float, default=None, help="Maximum inner-planet semimajor axis in AU. Defaults to observed planet periastron.")
    parser.add_argument("--da", type=float, default=0.001)
    parser.add_argument("--emin", type=float, default=0.000)
    parser.add_argument("--emax", type=float, default=0.500)
    parser.add_argument("--de", type=float, default=0.010)
    parser.add_argument("--inner-mass-earth", type=float, default=None, help="If omitted, uses observed CBP mass.")
    parser.add_argument("--orbits", type=float, default=1e5)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--jobid", type=int, default=0)
    parser.add_argument("--njobs", type=int, default=1)
    parser.add_argument("--outdir", type=str, default="CBP_out")
    return parser.parse_args()


def main():
    args = parse_args()
    sysdata = KEPLER_SYSTEMS[args.system]
    ensure_dir(args.outdir)
    amin = args.amin if args.amin is not None else 1.01 * sysdata["a_bin"]
    amax = args.amax if args.amax is not None else sysdata["a_planet"] * (1.0 - sysdata["e_planet"])
    inner_mass = normalize_mass_to_binary(args.inner_mass_earth * M_E if args.inner_mass_earth is not None else sysdata["m_planet"], sysdata)
    a_values = np.round(np.arange(amin, amax + 0.5 * args.da, args.da), 3)
    e_values = np.round(np.arange(args.emin, args.emax + 0.5 * args.de, args.de), 3)
    parameters = shuffled_parameters([(float(a), float(e)) for a in a_values for e in e_values])
    parameters = split_job(parameters, args.jobid, args.njobs)
    outer = {"a": sysdata["a_planet"] / sysdata["a_bin"], "e": sysdata["e_planet"], "m": normalize_mass_to_binary(sysdata["m_planet"], sysdata), "M": 180.0, "inc": 0.0}
    tag = args.system.replace("-", "")
    fname = f"{args.outdir}/MaxEcc_2pl_{tag}_job{args.jobid:02d}.txt" if args.njobs > 1 else f"{args.outdir}/MaxEcc_2pl_{tag}.txt"
    write_header(fname, "#a_0,e_0,emax,emin,tcol")
    def worker(par):
        a_au, ep = par
        a_ratio = a_au / sysdata["a_bin"]
        result = evaluate_sim(sysdata["mu"], sysdata["eb"], a_ratio, 0.0, ep=ep, inc=0.0, planet_mass=inner_mass, tscale=args.orbits, nout=300, a_bin=1.0, outer_planet=outer)
        result = result.__class__(a_au, result.m0, ep, result.emax, result.emin, result.tcol, result.status)
        return result_to_maxecc_row(result, include_e0=True)
    t0 = time.time()
    rows = run_pool(worker, parameters, workers=args.workers)
    append_rows(fname, rows)
    print(f"completed {len(rows)} integrations in {(time.time() - t0) / 60.0:.3f} min")
    print(f"wrote {fname}")


if __name__ == "__main__":
    main()
