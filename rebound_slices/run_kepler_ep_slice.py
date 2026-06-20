#!/usr/bin/env python3
"""Known-CBP host slice over test-planet semimajor axis and eccentricity.

Output follows the original CBP pool style:
#a_0,e_0,emax,emin,tcol
"""

import argparse
import time

import numpy as np

from cbp_common import KEPLER_SYSTEMS, M_J, append_rows, ensure_dir, evaluate_sim, normalize_mass_to_binary, result_to_maxecc_row, run_pool, shuffled_parameters, split_job, write_header


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", type=str, default="Kepler-16b", choices=sorted(KEPLER_SYSTEMS))
    parser.add_argument("--amin", type=float, default=None, help="Minimum test planet semimajor axis in AU.")
    parser.add_argument("--amax", type=float, default=None, help="Maximum test planet semimajor axis in AU.")
    parser.add_argument("--da", type=float, default=0.010)
    parser.add_argument("--emin", type=float, default=0.000)
    parser.add_argument("--emax", type=float, default=0.500)
    parser.add_argument("--de", type=float, default=0.010)
    parser.add_argument("--mass-jup", type=float, default=1.000)
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
    amax = args.amax if args.amax is not None else 1.50 * sysdata["a_planet"]
    a_values = np.round(np.arange(amin, amax + 0.5 * args.da, args.da), 3)
    e_values = np.round(np.arange(args.emin, args.emax + 0.5 * args.de, args.de), 3)
    parameters = shuffled_parameters([(float(a), float(e)) for a in a_values for e in e_values])
    parameters = split_job(parameters, args.jobid, args.njobs)
    tag = args.system.replace("-", "")
    fname = f"{args.outdir}/MaxEcc_{tag}_job{args.jobid:02d}.txt" if args.njobs > 1 else f"{args.outdir}/MaxEcc_{tag}.txt"
    write_header(fname, "#a_0,e_0,emax,emin,tcol")
    def worker(par):
        a_au, ep = par
        a_ratio = a_au / sysdata["a_bin"]
        result = evaluate_sim(sysdata["mu"], sysdata["eb"], a_ratio, 0.0, ep=ep, inc=0.0, planet_mass=normalize_mass_to_binary(args.mass_jup * M_J, sysdata), tscale=args.orbits, nout=300, a_bin=1.0)
        result = result.__class__(a_au, result.m0, ep, result.emax, result.emin, result.tcol, result.status)
        return result_to_maxecc_row(result, include_e0=True)
    t0 = time.time()
    rows = run_pool(worker, parameters, workers=args.workers)
    append_rows(fname, rows)
    print(f"completed {len(rows)} integrations in {(time.time() - t0) / 60.0:.3f} min")
    print(f"wrote {fname}")


if __name__ == "__main__":
    main()
