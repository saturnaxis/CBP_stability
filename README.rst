Stability Limit for the Innermost Circumbinary Planet
--------
.. image:: https://img.shields.io/badge/arXiv-1802.08868-red.svg?style=flat 
    :target: https://arxiv.org/abs/1802.08868 
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1174228.svg
   :target: https://doi.org/10.5281/zenodo.1174228
--------

This is a repository for tools to determine the stability limit, or smallest semimajor axis ratio (a_c = a_p/a_bin), of a circumbinary planet.  This summarized data in a_crit.txt collates ~150 million full N-body simulations using the well-established Mercury6 code (Chambers 1999, Chambers et al. 2002).  The simulations evaluated a wide range of initial binary parameters, including the binary mass ratio, mu, and binary eccentricity, e_bin.  The parameter mu is defined as the ratio of the smaller mass star divided by the total mass of the stellar binary.  These simulations were performed beginning at binary periastron and should be considered as conservative estimates of stability.  Due to symmetry only phases between 0 -- 180 degrees are tried, where the negative longitudes (-180 up to 0 deg.) are assumed to be mirror image of the positive ones.

After cloning the repository, the tool for determining a_c for a given set of binary parameters (mu, e_bin) can be determined simply by running **'python get_ac.py mu eb'**, where mu is replaced with a float [0, 0.5] and eb is replaced with a float [0,0.8].

The results from all the simulations are available on Zenodo.org (http://doi.org/10.5281/zenodo.1174228) as a compressed tar archive.  After downloading this tar archive, simply running **'python plot_from_tar.py mu eb'** produces a figure of the the given (mu, e_bin) parameters.  The value for mu must be chosen from 0.01 - 0.50 in steps of 0.01 and the value for e_bin must be chosen from 0.0 - 0.8 in steps of 0.01.

Attribution
--------
A more detailed description of these simulations, the use of an interpolation method, and the context for the Kepler circumbinary planets (CBPs) will be available in a forthcoming paper.  However in the meantime please use the following citation, if you find these tools useful in your research. ::

  @article{Quarles2018,
  author = {{Quarles}, B. and {Satyal}, S. and {Kostov}, V. and {Kaib}, N. and {Haghighipour}, N.},
  title = "{Stability Limits of Circumbinary Planets: Is There a Pile-up in the Kepler CBPs?}",
  journal = {ArXiv e-prints},
  archivePrefix = "arXiv",
  eprint = {1802.08868},
  primaryClass = "astro-ph.EP",
  keywords = {Astrophysics - Earth and Planetary Astrophysics},
  year = 2018,
  month = feb,
  adsurl = {http://adsabs.harvard.edu/abs/2018arXiv180208868Q},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
  }
