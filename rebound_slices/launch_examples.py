#!/usr/bin/env python3
"""Print example commands for local or SLURM-style slice launches."""

import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["hw", "inclination", "kepler", "packing", "beta"], default="hw")
    parser.add_argument("--njobs", type=int, default=10)
    return parser.parse_args()


def main():
    args = parse_args()
    for jobid in range(args.njobs):
        if args.mode == "hw":
            print(f"python run_hw_slice.py --mu 0.300 --eb 0.200 --jobid {jobid} --njobs {args.njobs} --workers 4")
        elif args.mode == "inclination":
            print(f"python run_inclination_slice.py --mu 0.300 --eb 0.200 --inc 45 --ac 2.800 --jobid {jobid} --njobs {args.njobs} --workers 4")
        elif args.mode == "kepler":
            print(f"python run_kepler_ep_slice.py --system Kepler-16b --jobid {jobid} --njobs {args.njobs} --workers 4")
        elif args.mode == "packing":
            print(f"python run_two_planet_packing_slice.py --system Kepler-16b --jobid {jobid} --njobs {args.njobs} --workers 4")
        elif args.mode == "beta":
            print(f"python run_beta_packing_slice.py --system Kepler-16b --jobid {jobid} --njobs {args.njobs} --workers 4")


if __name__ == "__main__":
    main()
