#!/usr/bin/env python3
"""Plot a beta-spacing two-planet slice from run_beta_packing_slice.py."""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cbp_common import ensure_dir


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str)
    parser.add_argument("--orbits", type=float, default=1e5)
    parser.add_argument("--outdir", type=str, default="figures")
    return parser.parse_args()


def main():
    args = parse_args()
    ensure_dir(args.outdir)
    data = pd.read_csv(args.input, comment="#", names=["beta", "e0", "emax", "emin", "tcol"])
    stable = np.isclose(data["tcol"], args.orbits)
    pivot = data.pivot_table(index="e0", columns="beta", values="emax", aggfunc="mean")
    stab = data.assign(stable=stable).pivot_table(index="e0", columns="beta", values="stable", aggfunc="mean")
    z = pivot.values.astype(float)
    z[stab.values < 1.0] = np.nan
    z = np.log10(np.maximum(z, 1e-3))
    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=180)
    im = ax.imshow(z, origin="lower", aspect="auto", extent=[pivot.columns.min(), pivot.columns.max(), pivot.index.min(), pivot.index.max()], vmin=-3, vmax=0)
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$e_0$")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\log_{10}(e_{max})$")
    base = os.path.splitext(os.path.basename(args.input))[0]
    outfile = os.path.join(args.outdir, base + "_emap.png")
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)
    print(f"wrote {outfile}")


if __name__ == "__main__":
    main()
