#!/usr/bin/env python3
"""Run one Quarles+2018-style single-planet stability slice using REBOUND.

Original-style output:
# a_0,M_0,emax,emin,tcol

Here a_0 is a_p/a_bin, M_0 is the initial planetary mean anomaly in degrees,
emax/emin are rounded to three decimals, and tcol is the instability time in
binary orbits. A negative tcol indicates an ejection-style escape.
"""

import argparse
import time

import numpy as np

from cbp_common import M_J, append_rows, ensure_dir, evaluate_sim, result_to_maxecc_row, run_pool, shuffled_parameters, split_job, write_header


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mu", type=float, default=0.300, help="Binary mass ratio mB/(mA+mB).")
    parser.add_argument("--eb", type=float, default=0.200, help="Binary eccentricity.")
    parser.add_argument("--amin", type=float, default=1.010, help="Minimum a_p/a_bin.")
    parser.add_argument("--amax", type=float, default=5.000, help="Maximum a_p/a_bin.")
    parser.add_argument("--da", type=float, default=0.010, help="Step in a_p/a_bin.")
    parser.add_argument("--mmin", type=float, default=0.000, help="Minimum initial mean anomaly in degrees.")
    parser.add_argument("--mmax", type=float, default=180.000, help="Maximum initial mean anomaly in degrees.")
    parser.add_argument("--dm", type=float, default=2.000, help="Step in mean anomaly in degrees.")
    parser.add_argument("--orbits", type=float, default=1e5, help="Integration time in binary orbits.")
    parser.add_argument("--nout", type=int, default=250, help="Number of output/check times.")
    parser.add_argument("--workers", type=int, default=1, help="Multiprocessing workers.")
    parser.add_argument("--jobid", type=int, default=0, help="Job slice index.")
    parser.add_argument("--njobs", type=int, default=1, help="Total number of job slices.")
    parser.add_argument("--seed", type=int, default=21417, help="Shuffle seed, matching the original style.")
    parser.add_argument("--outdir", type=str, default="GenRuns_out", help="Output folder.")
    parser.add_argument("--binary-phase", type=float, default=0.0, help="Initial binary true anomaly: 0=periastron, 180=apastron.")
    return parser.parse_args()


def main():
    args = parse_args()
    ensure_dir(args.outdir)
    a_values = np.round(np.arange(args.amin, args.amax + 0.5 * args.da, args.da), 3)
    m_values = np.round(np.arange(args.mmin, args.mmax + 0.5 * args.dm, args.dm), 3)
    parameters = [(float(a), float(m)) for a in a_values for m in m_values]
    parameters = shuffled_parameters(parameters, args.seed)
    parameters = split_job(parameters, args.jobid, args.njobs)
    fname = f"{args.outdir}/MaxEcc_[{args.mu:.3f},{args.eb:.3f}]_job{args.jobid:02d}.txt" if args.njobs > 1 else f"{args.outdir}/MaxEcc_[{args.mu:.3f},{args.eb:.3f}].txt"
    write_header(fname, "#a_0,M_0,emax,emin,tcol")
    def worker(par):
        a0, m0 = par
        result = evaluate_sim(args.mu, args.eb, a0, m0, ep=0.0, inc=0.0, planet_mass=M_J, tscale=args.orbits, nout=args.nout, binary_phase_deg=args.binary_phase)
        return result_to_maxecc_row(result, include_e0=False)
    t0 = time.time()
    rows = run_pool(worker, parameters, workers=args.workers)
    append_rows(fname, rows)
    print(f"completed {len(rows)} integrations in {(time.time() - t0) / 60.0:.3f} min")
    print(f"wrote {fname}")


if __name__ == "__main__":
    main()
