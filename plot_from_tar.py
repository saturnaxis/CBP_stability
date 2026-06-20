#!/usr/bin/env python3
"""Plot one stability map directly from the Zenodo ``MaxEcc.tar.gz`` archive."""

from __future__ import annotations

import argparse
from pathlib import Path

from cbp_plotting import plot_archive_two_panel, plot_phase_map, read_maxecc_from_tar


def main() -> None:
    parser = argparse.ArgumentParser(description='Plot a Quarles et al. (2018) MaxEcc grid from MaxEcc.tar.gz.')
    parser.add_argument('mu', type=float, help='Binary mass ratio on the archive grid, e.g. 0.200.')
    parser.add_argument('e_bin', type=float, help='Binary eccentricity on the archive grid, e.g. 0.200.')
    parser.add_argument('--tar', default='MaxEcc.tar.gz', help='Path to the Zenodo archive. Default: MaxEcc.tar.gz')
    parser.add_argument('--mode', choices=['emax', 'time', 'delta_e', 'two-panel'], default='emax')
    parser.add_argument('--output', default=None, help='Output PNG filename. A default name is used if omitted.')
    parser.add_argument('--dpi', type=int, default=300)
    args = parser.parse_args()

    data = read_maxecc_from_tar(args.tar, args.mu, args.e_bin)
    if args.output is None:
        suffix = 'archive_two_panel' if args.mode == 'two-panel' else args.mode
        args.output = f'MaxEcc_[{args.mu:0.3f},{args.e_bin:0.3f}]_{suffix}.png'
    if args.mode == 'two-panel':
        plot_archive_two_panel(data, output=args.output, dpi=args.dpi)
    else:
        plot_phase_map(data, output=args.output, title=fr'$\mu={args.mu:.3f}$, $e_{{bin}}={args.e_bin:.3f}$', mode=args.mode, dpi=args.dpi)
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
