#!/usr/bin/env python3
"""Two-planet beta-spacing slice for a Kepler CBP system.

This follows the logic of the original Mercury_CBP_beta.py script: fix an outer
planet at ao, choose a mutual-Hill spacing beta, convert beta to an interior
semimajor axis, and vary the shared initial eccentricity. The default phasing is
inner planet M=180 deg and outer planet M=0 deg, matching the original style.

Output:
#beta,e_0,emax,emin,tcol
"""

import argparse
import time

import numpy as np

from cbp_common import KEPLER_SYSTEMS, M_J, append_rows, ensure_dir, evaluate_sim, normalize_mass_to_binary, run_pool, shuffled_parameters, split_job, write_header


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", type=str, default="Kepler-16b", choices=sorted(KEPLER_SYSTEMS))
    parser.add_argument("--ao", type=float, default=None, help="Outer planet semimajor axis in AU. Defaults to the observed CBP semimajor axis.")
    parser.add_argument("--beta-min", type=float, default=1.000)
    parser.add_argument("--beta-max", type=float, default=25.000)
    parser.add_argument("--dbeta", type=float, default=0.005)
    parser.add_argument("--emin", type=float, default=0.000)
    parser.add_argument("--emax", type=float, default=0.500)
    parser.add_argument("--de", type=float, default=0.010)
    parser.add_argument("--planet-mass-jup", type=float, default=1.000, help="Mass of both planets, in Jupiter masses.")
    parser.add_argument("--orbits", type=float, default=1e5)
    parser.add_argument("--nout", type=int, default=300)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--jobid", type=int, default=0)
    parser.add_argument("--njobs", type=int, default=1)
    parser.add_argument("--seed", type=int, default=21417)
    parser.add_argument("--outdir", type=str, default="CBP_out")
    return parser.parse_args()


def inner_a_from_beta(ao_ratio, beta, planet_mass_norm):
    x = 0.5 * ((2.0 * planet_mass_norm) / 3.0) ** (1.0 / 3.0)
    return ao_ratio * (1.0 - beta * x) / (1.0 + beta * x)


def main():
    args = parse_args()
    sysdata = KEPLER_SYSTEMS[args.system]
    ensure_dir(args.outdir)
    ao_au = args.ao if args.ao is not None else sysdata["a_planet"]
    ao_ratio = ao_au / sysdata["a_bin"]
    m_norm = normalize_mass_to_binary(args.planet_mass_jup * M_J, sysdata)
    beta_values = np.round(np.arange(args.beta_min, args.beta_max + 0.5 * args.dbeta, args.dbeta), 3)
    e_values = np.round(np.arange(args.emin, args.emax + 0.5 * args.de, args.de), 3)
    parameters = shuffled_parameters([(float(beta), float(e)) for beta in beta_values for e in e_values], args.seed)
    parameters = split_job(parameters, args.jobid, args.njobs)
    tag = args.system.replace("-", "")
    fname = f"{args.outdir}/Beta_{tag}_job{args.jobid:02d}.txt" if args.njobs > 1 else f"{args.outdir}/Beta_{tag}.txt"
    write_header(fname, "#beta,e_0,emax,emin,tcol")
    outer = {"a": ao_ratio, "e": 0.0, "m": m_norm, "M": 0.0, "inc": 0.0}

    def worker(par):
        beta, ep = par
        a_inner = inner_a_from_beta(ao_ratio, beta, m_norm)
        peri_inner = a_inner * (1.0 - ep)
        if a_inner <= 0.0 or peri_inner <= 1.0:
            return f"{beta:.3f},{ep:.3f},{-1.000:.3f},{-1.000:.3f},{0.000:.3f}\n"
        result = evaluate_sim(sysdata["mu"], sysdata["eb"], a_inner, 180.0, ep=ep, inc=0.0, planet_mass=m_norm, tscale=args.orbits, nout=args.nout, a_bin=1.0, outer_planet=outer)
        return f"{beta:.3f},{ep:.3f},{result.emax:.3f},{result.emin:.3f},{result.tcol:.3f}\n"

    t0 = time.time()
    rows = run_pool(worker, parameters, workers=args.workers)
    append_rows(fname, rows)
    print(f"completed {len(rows)} integrations in {(time.time() - t0) / 60.0:.3f} min")
    print(f"wrote {fname}")


if __name__ == "__main__":
    main()
