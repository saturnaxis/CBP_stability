# Example notebooks

Run these notebooks from the repository root. They are ordered from the lightest, interpolation-only workflow to the heavier Zenodo and REBOUND workflows.

| Notebook | Purpose | Data needed |
|---|---|---|
| `00_interpolation_quickstart.ipynb` | Query `a_c` and convert to a critical period. | `a_crit.txt` |
| `01_figures_1_2_stability_maps.ipynb` | Recreate Figure 1/2 style stability maps from included subsets. | `plot_figures/Fig1_data`, `plot_figures/Fig2_data` |
| `02_figures_3_4_acrit_trends.ipynb` | Recreate median/range stability-limit trends. | `a_crit.txt` |
| `03_figure_5_stability_surface.ipynb` | Plot the stability surface and Kepler CBP locations. | `a_crit.txt` |
| `04_figure_6_packing_summary.ipynb` | Recreate the dynamical-spacing summary. | Embedded table values |
| `05_figures_7_8_kepler_maps.ipynb` | Explore Kepler-system single-planet and two-planet maps. | `plot_figures/Fig8_data` |
| `06_zenodo_archive_custom_maps.ipynb` | Plot arbitrary archive maps from `MaxEcc.tar.gz`. | Zenodo archive |
| `07_rebound_slice_recreation.ipynb` | Rerun small slices using REBOUND. | Optional `rebound` install |
