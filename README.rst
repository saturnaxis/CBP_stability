Stability Limits of Circumbinary Planets
========================================

.. image:: https://img.shields.io/badge/arXiv-1802.08868-red.svg?style=flat
   :target: https://arxiv.org/abs/1802.08868
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1174228.svg
   :target: https://doi.org/10.5281/zenodo.1174228
.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/saturnaxis/CBP_stability/master?filepath=CBP_stability_interpolation_example.ipynb

This repository provides data and tools associated with Quarles et al. (2018),
*Stability Limits of Circumbinary Planets: Is There a Pile-up in the Kepler
CBPs?*  The main quantity is the critical circumbinary semimajor-axis ratio

::

   a_c = a_p/a_bin,

where ``a_c`` is the smallest initially circular, coplanar planetary orbit that
survived the grid of initial planetary phases for a given binary mass ratio
``mu`` and binary eccentricity ``e_bin``.

The summarized grid in ``a_crit.txt`` comes from about 150 million Mercury6
N-body simulations. The full simulation archive is hosted on Zenodo as
``MaxEcc.tar.gz``. This repository now has three complementary layers:

1. **Interpolation tools** for quick stability-limit lookup using ``a_crit.txt``.
2. **Plotting tools and notebooks** for recreating paper-style figures and
   plotting arbitrary maps from the Zenodo archive.
3. **Optional REBOUND slice scripts** for rerunning small, modern validation
   slices of the original experiment. These are not replacements for the
   Mercury6 production archive.

Installation
------------

A minimal local setup is::

   conda env create -f environment.yml
   conda activate cbp-stability

or, with pip::

   pip install numpy scipy pandas matplotlib jupyter

The optional REBOUND examples additionally require::

   pip install rebound

Quick stability lookup
----------------------

Use ``get_ac.py`` from the repository root::

   python get_ac.py 0.230 0.159

which prints, for example::

   a_c = 2.880

If the binary period is known, add ``--p-bin-days`` to estimate the critical
planetary period using Kepler's third law::

   python get_ac.py 0.230 0.159 --p-bin-days 41.07758

The mass ratio is defined as

::

   mu = M_B/(M_A + M_B),

where ``M_B`` is the lower-mass stellar component. If observations quote
``q = M_B/M_A``, convert using::

   mu = q/(1 + q)

The interpolation is implemented in ``cbp_stability.py`` using SciPy's modern
``RegularGridInterpolator`` rather than the deprecated ``interp2d`` interface.

Notebook guide
--------------

The notebook examples are the recommended entry point for new users.

+------------------------------------------------+-----------------------------------------------+------------------------------+
| Notebook                                       | Purpose                                       | Data needed                  |
+================================================+===============================================+==============================+
| ``CBP_stability_interpolation_example.ipynb``  | Binder-friendly interpolation quickstart.     | ``a_crit.txt``               |
+------------------------------------------------+-----------------------------------------------+------------------------------+
| ``examples/00_interpolation_quickstart.ipynb`` | Minimal ``a_c`` and ``P_c`` example.           | ``a_crit.txt``               |
+------------------------------------------------+-----------------------------------------------+------------------------------+
| ``examples/01_figures_1_2_stability_maps.ipynb`` | Figure 1/2 style phase maps.                | Included Fig1/Fig2 subsets   |
+------------------------------------------------+-----------------------------------------------+------------------------------+
| ``examples/02_figures_3_4_acrit_trends.ipynb`` | Median/range trends of ``a_c`` vs ``e_bin``.  | ``a_crit.txt``               |
+------------------------------------------------+-----------------------------------------------+------------------------------+
| ``examples/03_figure_5_stability_surface.ipynb`` | ``a_c(mu,e_bin)`` surface and Kepler CBPs.  | ``a_crit.txt``               |
+------------------------------------------------+-----------------------------------------------+------------------------------+
| ``examples/04_figure_6_packing_summary.ipynb`` | Dynamical-spacing / planet-packing summary.  | Embedded table values        |
+------------------------------------------------+-----------------------------------------------+------------------------------+
| ``examples/05_figures_7_8_kepler_maps.ipynb``  | Kepler-system single/two-planet maps.         | Included Fig8 data subsets   |
+------------------------------------------------+-----------------------------------------------+------------------------------+
| ``examples/06_zenodo_archive_custom_maps.ipynb`` | Plot arbitrary maps from ``MaxEcc.tar.gz``. | Zenodo archive               |
+------------------------------------------------+-----------------------------------------------+------------------------------+
| ``examples/07_rebound_slice_recreation.ipynb`` | Rerun small validation slices with REBOUND.   | Optional REBOUND install     |
+------------------------------------------------+-----------------------------------------------+------------------------------+

Using the Zenodo archive
------------------------

Download ``MaxEcc.tar.gz`` from Zenodo and place it in the repository root.
Then plot any grid point in the archive without extracting the full tarball::

   python plot_from_tar.py 0.200 0.200 --mode emax

The available archive grid uses ``mu = 0.01`` to ``0.50`` in steps of ``0.01``
plus the low-mass ``mu = 0.001`` case, and ``e_bin = 0.00`` to ``0.80`` in
steps of ``0.01``. The script expects exact grid values.

The ``--mode`` option controls what is plotted:

+----------------+---------------------------------------------------+
| Mode           | Meaning                                           |
+================+===================================================+
| ``emax``       | Map of ``log10(e_max)`` with unstable cells blank. |
+----------------+---------------------------------------------------+
| ``time``       | Map of ``log10(|t_col|)``.                         |
+----------------+---------------------------------------------------+
| ``delta_e``    | Map of ``log10(e_max-e_min)``.                     |
+----------------+---------------------------------------------------+
| ``two-panel``  | Later diagnostic: lifetime map plus ``Delta e``.   |
+----------------+---------------------------------------------------+

If you have already extracted a text file from the archive, use::

   python plot_from_file.py MaxEcc_[0.200,0.200].txt --mode two-panel

The original first-update plotting scripts are preserved in
``original_later_plotting/`` for reference. The modern top-level scripts use the
same ideas, but are Python 3 compatible and share helper functions in
``cbp_plotting.py``.

Included figure data
--------------------

The repository includes reduced data subsets for several paper figures:

``plot_figures/Fig1_data``
   Periastron-start 4-by-4 phase-map subset.

``plot_figures/Fig2_data``
   Apastron-start comparison subset.

``plot_figures/Fig8_data``
   Kepler-system single-planet and two-planet map subsets.

These files make the example notebooks useful without requiring the full Zenodo
archive.

Optional REBOUND slice scripts
------------------------------

The ``rebound_slices/`` directory contains small, portable scripts that recreate
pieces of the original experiment with REBOUND. These scripts follow the style
of the original project scripts: grid-based jobs, comma-separated output files,
job splitting, and separate reduction/plotting scripts. They are intended for
validation, teaching, and exploratory extensions.

A tiny smoke test is::

   python rebound_slices/run_hw_slice.py --mu 0.300 --eb 0.200 --amin 2.0 --amax 2.2 --da 0.1 --dm 30 --orbits 100 --workers 1
   python rebound_slices/reduce_hw_acrit.py --orbits 100
   python rebound_slices/plot_hw_slice.py GenRuns_out/MaxEcc_[0.300,0.200].txt --orbits 100

For production-like comparisons, use the paper's much longer integration time
of ``10^5`` binary orbits and the original grid spacing. Expect small boundary
differences because these examples use REBOUND rather than the original Mercury6
production setup.

Repository layout
-----------------

::

   a_crit.txt                         # Published stability-limit summary grid
   cbp_stability.py                   # Interpolation utilities
   cbp_plotting.py                    # Plotting utilities for MaxEcc files
   get_ac.py                          # Command-line stability lookup
   plot_from_tar.py                   # Plot maps directly from MaxEcc.tar.gz
   plot_from_file.py                  # Plot maps from extracted MaxEcc text files
   CBP_stability_interpolation_example.ipynb
   examples/                          # Notebook guide to common workflows
   plot_figures/                      # Included paper-figure data subsets
   original_later_plotting/           # Original post-paper plotting scripts
   rebound_slices/                    # Optional REBOUND recreation scripts

Citation
--------

If you use these tools or data, please cite::

  @article{Quarles2018,
    author = {{Quarles}, B. and {Satyal}, S. and {Kostov}, V. and {Kaib}, N. and {Haghighipour}, N.},
    title = "{Stability Limits of Circumbinary Planets: Is There a Pile-up in the Kepler CBPs?}",
    journal = {The Astrophysical Journal},
    year = 2018,
    month = apr,
    volume = 856,
    eid = {150},
    pages = {150},
    doi = {10.3847/1538-4357/aab264},
    archivePrefix = {arXiv},
    eprint = {1802.08868},
    primaryClass = {astro-ph.EP}
  }
