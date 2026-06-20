#!/usr/bin/env python3
"""Make a 4x4 paper-style panel from MaxEcc slice files.

This is the REBOUND analog of the original plot_muEcc.py / plot_muEcc_emax.py
scripts. It expects files named like MaxEcc_[0.300,0.200].txt, with columns
#a_0,M_0,emax,emin,tcol.
"""

import argparse
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cbp_common import ensure_dir


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", type=str, default="GenRuns_out")
    parser.add_argument("--outdir", type=str, default="figures")
    parser.add_argument("--mode", choices=["emax", "stable"], default="emax")
    parser.add_argument("--orbits", type=float, default=1e5)
    parser.add_argument("--mu-values", type=str, default="0.001,0.100,0.300,0.500")
    parser.add_argument("--eb-values", type=str, default="0.000,0.100,0.300,0.500")
    parser.add_argument("--outfile", type=str, default=None)
    return parser.parse_args()


def parse_mu_eb(path):
    match = re.search(r"MaxEcc_\[(?P<mu>[0-9.]+),(?P<eb>[0-9.]+)\]", os.path.basename(path))
    if match is None:
        return np.nan, np.nan
    return float(match.group("mu")), float(match.group("eb"))


def find_file(indir, mu, eb):
    pattern = os.path.join(indir, glob.escape(f"MaxEcc_[{mu:.3f},{eb:.3f}]") + "*.txt")
    files = sorted(f for f in glob.glob(pattern) if "Incl" not in os.path.basename(f) and "job" not in os.path.basename(f))
    if files:
        return files[0]
    all_files = glob.glob(os.path.join(indir, "MaxEcc_[[]*.txt"))
    for path in all_files:
        this_mu, this_eb = parse_mu_eb(path)
        if np.isclose(this_mu, mu, atol=5e-4) and np.isclose(this_eb, eb, atol=5e-4) and "Incl" not in os.path.basename(path):
            return path
    return None


def acrit_from_data(data, orbits):
    work = data.copy()
    work["stable"] = np.isclose(work["tcol"], orbits)
    grouped = work.groupby("a0")["stable"].all().reset_index(name="stable")
    stable_all = grouped[grouped["stable"]]
    if len(stable_all) == 0:
        return np.nan
    return float(stable_all["a0"].min())


def main():
    args = parse_args()
    ensure_dir(args.outdir)
    mu_values = [float(x) for x in args.mu_values.split(",")]
    eb_values = [float(x) for x in args.eb_values.split(",")]
    fig, axes = plt.subplots(len(mu_values), len(eb_values), figsize=(13, 13), dpi=180, sharex=True, sharey=True)
    if len(mu_values) == 1 and len(eb_values) == 1:
        axes = np.array([[axes]])
    elif len(mu_values) == 1:
        axes = axes[None, :]
    elif len(eb_values) == 1:
        axes = axes[:, None]
    last_im = None
    for i, mu in enumerate(mu_values):
        for j, eb in enumerate(eb_values):
            ax = axes[i, j]
            path = find_file(args.indir, mu, eb)
            if path is None:
                ax.text(0.5, 0.5, "missing", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"({mu:.3f}, {eb:.3f})", fontsize=10)
                continue
            data = pd.read_csv(path, comment="#", names=["a0", "M0", "emax", "emin", "tcol"])
            ac = acrit_from_data(data, args.orbits)
            pivot = data.pivot_table(index="a0", columns="M0", values="emax", aggfunc="mean")
            stable = data.pivot_table(index="a0", columns="M0", values="tcol", aggfunc="mean")
            if args.mode == "stable":
                z = np.where(np.isclose(stable.values, args.orbits), 1.0, 0.0)
                vmin, vmax = 0.0, 1.0
                color_label = "stable"
            else:
                z = pivot.values.astype(float)
                z[~np.isclose(stable.values, args.orbits)] = np.nan
                z = np.log10(np.maximum(z, 1e-3))
                vmin, vmax = -3.0, 0.0
                color_label = r"$\log_{10}(e_{max})$"
            last_im = ax.imshow(z, origin="lower", aspect="auto", extent=[pivot.columns.min(), pivot.columns.max(), pivot.index.min(), pivot.index.max()], vmin=vmin, vmax=vmax)
            if np.isfinite(ac):
                ax.axhline(ac, lw=2)
                ax.text(0.03, 0.05, f"$a_c$={ac:.3f}", fontsize=8, transform=ax.transAxes)
            ax.text(0.97, 0.90, f"({mu:.3f}, {eb:.2f})", ha="right", fontsize=9, transform=ax.transAxes)
    fig.subplots_adjust(hspace=0.0, wspace=0.0, right=0.88)
    fig.text(0.5, 0.04, "Mean anomaly (deg)", ha="center")
    fig.text(0.04, 0.5, r"$a_p/a_{bin}$", va="center", rotation="vertical")
    if last_im is not None:
        cax = fig.add_axes([0.90, 0.12, 0.015, 0.76])
        cbar = fig.colorbar(last_im, cax=cax)
        cbar.set_label(color_label)
    outfile = args.outfile or os.path.join(args.outdir, f"MuEcc_bin0_{args.mode}.png")
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outfile}")


if __name__ == "__main__":
    main()
