#!/usr/bin/env python3
"""Run a small survey over (mu, e_bin), one grid per pair.

This mirrors the original launch/splitting idea but keeps each grid as a local
REBOUND job. For production, use small job slices on a cluster.
"""

import argparse
import subprocess
import sys

import numpy as np

from cbp_common import shuffled_parameters, split_job


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mu-min", type=float, default=0.010)
    parser.add_argument("--mu-max", type=float, default=0.500)
    parser.add_argument("--dmu", type=float, default=0.010)
    parser.add_argument("--eb-min", type=float, default=0.000)
    parser.add_argument("--eb-max", type=float, default=0.800)
    parser.add_argument("--deb", type=float, default=0.010)
    parser.add_argument("--jobid", type=int, default=0)
    parser.add_argument("--njobs", type=int, default=1)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--orbits", type=float, default=1e5)
    parser.add_argument("--amin", type=float, default=1.010)
    parser.add_argument("--amax", type=float, default=5.000)
    parser.add_argument("--da", type=float, default=0.010)
    parser.add_argument("--dm", type=float, default=2.000)
    parser.add_argument("--seed", type=int, default=21417)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    mu_values = np.round(np.arange(args.mu_min, args.mu_max + 0.5 * args.dmu, args.dmu), 3)
    eb_values = np.round(np.arange(args.eb_min, args.eb_max + 0.5 * args.deb, args.deb), 3)
    parameters = [(float(mu), float(eb)) for mu in mu_values for eb in eb_values]
    parameters = shuffled_parameters(parameters, args.seed)
    parameters = split_job(parameters, args.jobid, args.njobs)
    for mu, eb in parameters:
        cmd = [sys.executable, "run_hw_slice.py", "--mu", f"{mu:.3f}", "--eb", f"{eb:.3f}", "--workers", str(args.workers), "--orbits", str(args.orbits), "--amin", str(args.amin), "--amax", str(args.amax), "--da", str(args.da), "--dm", str(args.dm)]
        print(" ".join(cmd))
        if not args.dry_run:
            subprocess.check_call(cmd)


if __name__ == "__main__":
    main()
