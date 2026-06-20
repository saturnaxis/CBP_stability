#!/usr/bin/env python3
"""Plot one stability map from an extracted ``MaxEcc`` text file."""

from __future__ import annotations

import argparse
from pathlib import Path

from cbp_plotting import plot_archive_two_panel, plot_phase_map, read_maxecc_file


def main() -> None:
    parser = argparse.ArgumentParser(description='Plot a Quarles et al. (2018) MaxEcc text grid.')
    parser.add_argument('file', help='Path to MaxEcc_[mu,e_bin].txt or similar.')
    parser.add_argument('--mode', choices=['emax', 'time', 'delta_e', 'two-panel'], default='two-panel')
    parser.add_argument('--output', default=None, help='Output PNG filename. A default name is used if omitted.')
    parser.add_argument('--dpi', type=int, default=300)
    args = parser.parse_args()

    data = read_maxecc_file(args.file)
    if args.output is None:
        stem = Path(args.file).stem.replace('[', '').replace(']', '').replace(',', '_')
        args.output = f'{stem}_{args.mode}.png'
    if args.mode == 'two-panel':
        plot_archive_two_panel(data, output=args.output, dpi=args.dpi)
    else:
        plot_phase_map(data, output=args.output, mode=args.mode, dpi=args.dpi)
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
