#!/usr/bin/env python3
"""Fit analytic stability-limit coefficients from an a_crit table.

This is a portable, Python-3 version of the original fit_emcee.py workflow. The
default fit is the HW99/Q2018-style polynomial

    a_c = C1 + C2 e + C3 e^2 + C4 mu + C5 mu e + C6 mu^2 + C7 (mu e)^2.

Use --run-emcee to optionally sample coefficient uncertainties with emcee and
corner, if those packages are installed.
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from cbp_common import ensure_dir


def model(theta, mu, eb):
    c1, c2, c3, c4, c5, c6, c7 = theta
    return c1 + c2 * eb + c3 * eb**2 + c4 * mu + c5 * mu * eb + c6 * mu**2 + c7 * (mu * eb)**2


def nll(theta, mu, eb, ac, sigma):
    residual = (ac - model(theta, mu, eb)) / sigma
    return 0.5 * np.sum(residual**2 + np.log(2.0 * np.pi * sigma**2))


def log_prior(theta):
    c1, c2, c3, c4, c5, c6, c7 = theta
    bounds = [(1.0, 2.5), (2.0, 6.5), (-4.0, 1.0), (0.0, 7.0), (-8.0, 4.0), (-12.0, 1.0), (-10.0, 8.0)]
    for value, (lo, hi) in zip(theta, bounds):
        if not lo < value < hi:
            return -np.inf
    return 0.0


def log_prob(theta, mu, eb, ac, sigma):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp - nll(theta, mu, eb, ac, sigma)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--acrit", type=str, default="a_crit_rebound.txt")
    parser.add_argument("--sigma", type=float, default=0.005)
    parser.add_argument("--outdir", type=str, default="fit_outputs")
    parser.add_argument("--run-emcee", action="store_true")
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--burn", type=int, default=1000)
    parser.add_argument("--walkers", type=int, default=64)
    parser.add_argument("--seed", type=int, default=21417)
    return parser.parse_args()


def main():
    args = parse_args()
    ensure_dir(args.outdir)
    df = pd.read_csv(args.acrit, comment="#", names=["mu", "eb", "ac"])
    df = df[np.isfinite(df["ac"])]
    mu = df["mu"].to_numpy(float)
    eb = df["eb"].to_numpy(float)
    ac = df["ac"].to_numpy(float)
    sigma = np.ones_like(ac) * args.sigma
    theta0 = np.array([1.48, 3.92, -1.41, 5.14, 0.33, -7.95, -4.89], dtype=float)
    result = minimize(nll, theta0, args=(mu, eb, ac, sigma), method="Nelder-Mead")
    coeff = result.x
    pred = model(coeff, mu, eb)
    chi2 = np.sum(((pred - ac) / sigma) ** 2)
    dof = max(1, len(ac) - len(coeff))
    outfile = os.path.join(args.outdir, "acrit_fit_coefficients.txt")
    with open(outfile, "w", encoding="utf-8") as f:
        f.write("#C1,C2,C3,C4,C5,C6,C7,chi2,reduced_chi2\n")
        f.write(",".join(f"{x:.3f}" for x in coeff) + f",{chi2:.3f},{chi2 / dof:.3f}\n")
    print(f"wrote {outfile}")
    print("coefficients:", ", ".join(f"{x:.3f}" for x in coeff))
    print(f"reduced chi2: {chi2 / dof:.3f}")

    if args.run_emcee:
        try:
            import emcee
            import corner
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise SystemExit("--run-emcee requires emcee and corner") from exc
        rng = np.random.default_rng(args.seed)
        ndim = len(coeff)
        pos = coeff + 1e-4 * rng.normal(size=(args.walkers, ndim))
        sampler = emcee.EnsembleSampler(args.walkers, ndim, log_prob, args=(mu, eb, ac, sigma))
        sampler.run_mcmc(pos, args.burn, progress=False)
        sampler.reset()
        sampler.run_mcmc(None, args.steps, progress=False)
        samples = sampler.get_chain(flat=True)
        chain_path = os.path.join(args.outdir, "chain.dat")
        np.savetxt(chain_path, samples, fmt="%.8e")
        fig = corner.corner(samples, labels=[f"C{i}" for i in range(1, 8)], quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".3f")
        fig.savefig(os.path.join(args.outdir, "triangle.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {chain_path}")


if __name__ == "__main__":
    main()
