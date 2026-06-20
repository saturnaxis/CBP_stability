#!/usr/bin/env python3
"""Reduce MaxEcc files into an a_crit table using the paper definition.

For each (mu, e_bin) slice, a_crit is the smallest a_p/a_bin where every tested
initial planetary phase survives to the requested integration time. Job-split
files are grouped and concatenated before the stability limit is measured.
"""

import argparse
import glob
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd

from cbp_common import ensure_dir


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", type=str, default="GenRuns_out")
    parser.add_argument("--outfile", type=str, default="a_crit_rebound.txt")
    parser.add_argument("--orbits", type=float, default=1e5)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--include-inclination", action="store_true", help="Include files with _Incl[...] in their names.")
    return parser.parse_args()


def parse_mu_eb(path):
    match = re.search(r"MaxEcc_\[(?P<mu>[0-9.]+),(?P<eb>[0-9.]+)\]", os.path.basename(path))
    if match is None:
        return np.nan, np.nan
    return float(match.group("mu")), float(match.group("eb"))


def acrit_from_data(data, orbits, tol):
    if len(data) == 0:
        return np.nan
    data = data.drop_duplicates(subset=["a0", "M0"], keep="last").copy()
    data["stable"] = np.abs(data["tcol"] - orbits) <= tol
    grouped = data.groupby("a0")["stable"].all().reset_index()
    stable_all = grouped[grouped["stable"]]
    if len(stable_all) == 0:
        return np.nan
    return float(stable_all["a0"].min())


def main():
    args = parse_args()
    ensure_dir(os.path.dirname(args.outfile) or ".")
    files = sorted(glob.glob(os.path.join(args.indir, "MaxEcc_[[]*.txt")))
    groups = defaultdict(list)
    for path in files:
        if (not args.include_inclination) and "_Incl[" in os.path.basename(path):
            continue
        mu, eb = parse_mu_eb(path)
        if np.isfinite(mu):
            groups[(mu, eb)].append(path)
    rows = []
    for (mu, eb), paths in sorted(groups.items()):
        frames = []
        for path in paths:
            frames.append(pd.read_csv(path, comment="#", names=["a0", "M0", "emax", "emin", "tcol"]))
        data = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["a0", "M0", "emax", "emin", "tcol"])
        rows.append((mu, eb, acrit_from_data(data, args.orbits, args.tol)))
    with open(args.outfile, "w", encoding="utf-8") as f:
        f.write("#mu,eb,a_crit\n")
        for mu, eb, ac in rows:
            f.write(f"{mu:.3f},{eb:.3f},{ac:.3f}\n")
    print(f"wrote {args.outfile}")


if __name__ == "__main__":
    main()
