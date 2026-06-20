#!/usr/bin/env python3
"""Plot a REBOUND MaxEcc slice in the style of the Quarles+2018 maps."""

import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cbp_common import ensure_dir


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str, help="MaxEcc_[mu,eb].txt file.")
    parser.add_argument("--orbits", type=float, default=1e5)
    parser.add_argument("--outdir", type=str, default="figures")
    parser.add_argument("--show", action="store_true")
    return parser.parse_args()


def parse_mu_eb(path):
    match = re.search(r"MaxEcc_\[(?P<mu>[0-9.]+),(?P<eb>[0-9.]+)\]", os.path.basename(path))
    if match is None:
        return np.nan, np.nan
    return float(match.group("mu")), float(match.group("eb"))


def main():
    args = parse_args()
    ensure_dir(args.outdir)
    mu, eb = parse_mu_eb(args.input)
    data = pd.read_csv(args.input, comment="#", names=["a0", "M0", "emax", "emin", "tcol"])
    data["stable"] = np.isclose(data["tcol"], args.orbits)
    grouped = data.groupby("a0")["stable"].all().reset_index()
    stable_all = grouped[grouped["stable"]]
    ac = np.nan if len(stable_all) == 0 else float(stable_all["a0"].min())
    pivot = data.pivot_table(index="a0", columns="M0", values="emax", aggfunc="mean")
    stable = data.pivot_table(index="a0", columns="M0", values="stable", aggfunc="mean")
    z = pivot.values.astype(float)
    z[stable.values < 1.0] = np.nan
    z = np.log10(np.maximum(z, 1e-3))
    fig, ax = plt.subplots(figsize=(7.0, 5.0), dpi=180)
    im = ax.imshow(z, origin="lower", aspect="auto", extent=[pivot.columns.min(), pivot.columns.max(), pivot.index.min(), pivot.index.max()], vmin=-3, vmax=0)
    if np.isfinite(ac):
        ax.axhline(ac, lw=2)
        ax.text(0.02, 0.04, f"$a_c$ = {ac:.3f}", transform=ax.transAxes, fontsize=10)
    ax.set_xlabel("Mean anomaly (deg)")
    ax.set_ylabel(r"$a_p/a_{bin}$")
    ax.set_title(rf"$\mu={mu:.3f}$, $e_{{bin}}={eb:.3f}$")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\log_{10}(e_{max})$")
    outfile = os.path.join(args.outdir, f"MuEcc_[{mu:.3f},{eb:.3f}]_emax.png")
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    if args.show:
        plt.show()
    plt.close(fig)
    print(f"wrote {outfile}")


if __name__ == "__main__":
    main()
