Stability Limit for the Innermost Circumbinary Planet
--------

This is a repository for tools to determine the stability limit, or smallest semimajor axis ratio (a_c = a_p/a_bin), of a circumbinary planet.  This summarized data in a_crit.txt collates ~150 million full N-body simulations using the well-established Mercury6 code (Chambers 1999, Chambers et al. 2002).  The simulations evaluated a wide range of initial binary parameters, including the binary mass ratio, mu, and binary eccentricity, e_bin.  The parameter mu is defined as the ratio of the smaller mass star divided by the total mass of the stellar binary.  These simulations were performed beginning at binary periastron and should be considered as conservative estimates of stability.  Due to symmetry only phases between 0 -- 180 degrees are tried, where the negative longitudes (-180 up to 0 deg.) are assumed to be mirror image of the positive ones.

After cloning the repository, the tool for determining a_c for a given set of binary parameters (mu, e_bin) can be determined simply by running 'python get_ac.py mu eb', where the mu is replaced with a float [0, 0.5] and eb is replaced with a float [0,0.8].

Attribution
--------
A more detailed description of these simulations, the use of an interpolation method, and the context for the Kepler circumbinary planets (CBPs) will be available in a forthcoming paper.  However in the meantime please use the following citation, if you find these tools useful in your research.::

  @article{Quarles2018,
  author = {{Quarles}, B. and {Saytal}, S. and {Kostov}, V. and {Kaib}, N.} and {Haghighipour}, N.,
  title = {Stability Limits of Circumbinary Planets: Is There a Pile-up in the Kepler CBPs?},
  journal = {ApJ},
  year = 2018,
  pages = {under review}
  }
