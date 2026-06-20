#!/usr/bin/env python3
"""Command-line lookup for the Quarles et al. (2018) stability limit."""

from __future__ import annotations

import argparse

from cbp_stability import build_stability_grid, critical_period_days, get_ac


def main() -> None:
    parser = argparse.ArgumentParser(description='Interpolate the CBP stability limit a_c = a_p/a_bin.')
    parser.add_argument('mu', type=float, help='Binary mass ratio: smaller stellar mass divided by total stellar mass.')
    parser.add_argument('e_bin', type=float, help='Binary eccentricity.')
    parser.add_argument('--p-bin-days', type=float, default=None, help='Optional binary period in days. If supplied, print P_c too.')
    parser.add_argument('--decimals', type=int, default=3, help='Number of decimals to print. Default: 3.')
    args = parser.parse_args()

    grid = build_stability_grid()
    ac = get_ac(args.mu, args.e_bin, grid=grid)
    print(f'a_c = {ac:.{args.decimals}f}')
    if args.p_bin_days is not None:
        pc = critical_period_days(args.mu, args.e_bin, args.p_bin_days, grid=grid)
        print(f'P_c = {pc:.{args.decimals}f} days')


if __name__ == '__main__':
    main()
