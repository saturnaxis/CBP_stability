#!/usr/bin/env python3
"""Inclined single-planet slice, modeled after Gen_ICs_HW_ext.py."""

import argparse
import time

import numpy as np

from cbp_common import M_J, append_rows, ensure_dir, evaluate_sim, result_to_maxecc_row, run_pool, shuffled_parameters, split_job, write_header


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mu", type=float, default=0.300)
    parser.add_argument("--eb", type=float, default=0.200)
    parser.add_argument("--inc", type=float, default=45.000)
    parser.add_argument("--ac", type=float, default=None, help="Optional center for grid. If supplied, defaults to 0.5ac to 1.5ac.")
    parser.add_argument("--amin", type=float, default=None)
    parser.add_argument("--amax", type=float, default=None)
    parser.add_argument("--da", type=float, default=0.010)
    parser.add_argument("--dm", type=float, default=2.000)
    parser.add_argument("--orbits", type=float, default=5e5)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--jobid", type=int, default=0)
    parser.add_argument("--njobs", type=int, default=1)
    parser.add_argument("--outdir", type=str, default="GenRuns_out")
    return parser.parse_args()


def main():
    args = parse_args()
    ensure_dir(args.outdir)
    amin = args.amin if args.amin is not None else (0.5 * args.ac if args.ac is not None else 1.010)
    amax = args.amax if args.amax is not None else (1.5 * args.ac if args.ac is not None else 5.000)
    a_values = np.round(np.arange(amin, amax + 0.5 * args.da, args.da), 3)
    m_values = np.round(np.arange(0.0, 180.0 + 0.5 * args.dm, args.dm), 3)
    parameters = shuffled_parameters([(float(a), float(m)) for a in a_values for m in m_values])
    parameters = split_job(parameters, args.jobid, args.njobs)
    fname = f"{args.outdir}/MaxEcc_[{args.mu:.3f},{args.eb:.3f}]_Incl[{int(round(args.inc))}]_job{args.jobid:02d}.txt" if args.njobs > 1 else f"{args.outdir}/MaxEcc_[{args.mu:.3f},{args.eb:.3f}]_Incl[{int(round(args.inc))}].txt"
    write_header(fname, "#a_0,M_0,emax,emin,tcol")
    def worker(par):
        a0, m0 = par
        result = evaluate_sim(args.mu, args.eb, a0, m0, ep=0.0, inc=args.inc, planet_mass=M_J, tscale=args.orbits, nout=300)
        return result_to_maxecc_row(result, include_e0=False)
    t0 = time.time()
    rows = run_pool(worker, parameters, workers=args.workers)
    append_rows(fname, rows)
    print(f"completed {len(rows)} integrations in {(time.time() - t0) / 60.0:.3f} min")
    print(f"wrote {fname}")


if __name__ == "__main__":
    main()
