#!/usr/bin/env python3
"""Merge job-split MaxEcc/Beta outputs into one original-style file.

The run scripts write job-specific files when --njobs > 1. This utility removes
extra headers, concatenates rows, and sorts by the leading columns when possible.
"""

import argparse
import glob
import os
import re

import pandas as pd

from cbp_common import ensure_dir


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pattern", type=str, help="Glob pattern, e.g. 'GenRuns_out/MaxEcc_[0.300,0.200]_job*.txt'.")
    parser.add_argument("--outfile", type=str, default=None)
    return parser.parse_args()


def infer_outfile(files):
    first = files[0]
    base = os.path.basename(first)
    base = re.sub(r"_job\d+", "", base)
    return os.path.join(os.path.dirname(first), base)


def main():
    args = parse_args()
    files = sorted(glob.glob(args.pattern))
    if not files:
        bracket_safe = glob.escape(args.pattern).replace("[*]", "*").replace("[?]", "?")
        files = sorted(glob.glob(bracket_safe))
    if not files:
        raise SystemExit(f"no files matched {args.pattern}")
    outfile = args.outfile or infer_outfile(files)
    ensure_dir(os.path.dirname(outfile) or ".")
    header = None
    frames = []
    for path in files:
        with open(path, "r", encoding="utf-8") as f:
            first = f.readline().strip()
        if first.startswith("#") and header is None:
            header = first
        names = header.lstrip("#").split(",") if header is not None else None
        frames.append(pd.read_csv(path, comment="#", names=names))
    data = pd.concat(frames, ignore_index=True)
    sort_cols = [c for c in data.columns[:2]]
    data = data.sort_values(sort_cols).drop_duplicates(subset=sort_cols, keep="last")
    with open(outfile, "w", encoding="utf-8") as f:
        if header is not None:
            f.write(header + "\n")
        data.to_csv(f, index=False, header=False, float_format="%.3f")
    print(f"wrote {outfile}")


if __name__ == "__main__":
    main()
