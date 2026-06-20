#!/usr/bin/env python3
"""Plot the a_crit surface and a_c(e_bin) summary from an a_crit table.

This combines the roles of the original plot_acr.py and plot_acritEcc.py for the
REBOUND slice package. Input is a table with columns #mu,eb,a_crit.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cbp_common import KEPLER_SYSTEMS, ensure_dir


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--acrit", type=str, default="a_crit_rebound.txt")
    parser.add_argument("--outdir", type=str, default="figures")
    parser.add_argument("--prefix", type=str, default="rebound")
    parser.add_argument("--plot-kepler", action="store_true")
    return parser.parse_args()


def hw99_population_fit(eb):
    return 2.278 + 3.824 * eb - 1.71 * eb**2


def dvorak_fit(eb):
    return 2.37 + 2.76 * eb - 1.04 * eb**2


def main():
    args = parse_args()
    ensure_dir(args.outdir)
    df = pd.read_csv(args.acrit, comment="#", names=["mu", "eb", "ac"])
    df = df[np.isfinite(df["ac"])]
    df = df.sort_values(["mu", "eb"])

    pivot = df.pivot_table(index="eb", columns="mu", values="ac", aggfunc="mean")
    fig, ax = plt.subplots(figsize=(8.0, 5.4), dpi=180)
    im = ax.imshow(pivot.values, origin="lower", aspect="auto", extent=[pivot.columns.min(), pivot.columns.max(), pivot.index.min(), pivot.index.max()], vmin=np.nanmin(df["ac"]), vmax=np.nanmax(df["ac"]), interpolation="nearest")
    if args.plot_kepler:
        for name, sysdata in KEPLER_SYSTEMS.items():
            ax.plot(sysdata["mu"], sysdata["eb"], ".", color="white", ms=7)
            ax.text(sysdata["mu"], sysdata["eb"] + 0.015, name.replace("Kepler-", "").replace("b", ""), color="white", fontsize=7, ha="center")
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$e_{bin}$")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$a_c/a_{bin}$")
    surface_out = os.path.join(args.outdir, f"{args.prefix}_acrit_surface.png")
    fig.tight_layout()
    fig.savefig(surface_out, dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.4), dpi=180)
    mus = sorted(df["mu"].unique())
    for mu in mus:
        sub = df[np.isclose(df["mu"], mu)]
        if len(sub) > 0:
            ax.plot(sub["eb"], sub["ac"], ".", ms=3)
    stats = df[df["mu"] >= 0.1].groupby("eb")["ac"].agg(["median", "min", "max"]).reset_index()
    if len(stats) > 0:
        yerr = np.vstack([stats["median"] - stats["min"], stats["max"] - stats["median"]])
        ax.errorbar(stats["eb"], stats["median"], yerr=yerr, fmt="o", color="black", ms=4, label="median and range")
    x = np.linspace(max(0.0, df["eb"].min()), min(0.8, df["eb"].max()), 300)
    ax.plot(x, dvorak_fit(x), color="black", lw=1.5, label="Dvorak et al. style")
    ax.plot(x, hw99_population_fit(x), color="black", lw=1.5, ls="--", label="HW99 style")
    for n in range(2, 12):
        ax.axhline(n ** (2.0 / 3.0), color="0.75", lw=0.8, zorder=0)
    ax.set_xlabel(r"$e_{bin}$")
    ax.set_ylabel(r"$a_c/a_{bin}$")
    ax.set_xlim(-0.005, 0.805)
    ax.set_ylim(max(1.0, np.nanmin(df["ac"]) - 0.2), np.nanmax(df["ac"]) + 0.2)
    ax.legend(fontsize=8)
    eb_out = os.path.join(args.outdir, f"{args.prefix}_acrit_vs_eb.png")
    fig.tight_layout()
    fig.savefig(eb_out, dpi=300)
    plt.close(fig)
    print(f"wrote {surface_out}")
    print(f"wrote {eb_out}")


if __name__ == "__main__":
    main()
