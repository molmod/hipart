Using HiPart
============

For the terminology used in this chapter we refer to the :ref:`theory`. 

Preparing the atomic database
-----------------------------

Both the :ref:`hirshfeld` and :ref:`hirshfeld-i` rely on a database of
pro-atoms. These are spherically averaged atomic densities obtained from
single-atom computations. The program ``hi-atomdb.py`` generates such a dabase
using the Gaussian03 program.

The online help is as follows::

    toon@poony# hi-atomdb.py --help
    Usage: hi-atomdb.py [options] lot atoms

    hi-atomdb.py computes a database of pro-atomic densities.

    The following arguments are mandatory:
      * lot  --  The level of theory to be used in Gaussian03 input notation.
      * atoms  -- The atoms to be computed. One can specify ranges, e.g 1,2-5'
                  (avoid whitespace)

    It is recommended to run this script in a directory that is initially empty and
    that will contain nothing but the generated atom database. This script will
    generate quite a few files and subdirectories.

    Examples:

    hi-atomdb.py MP2/Aug-CC-pVDZ 1-10,17
    hi-atomdb.py HF/3-21G 1,6,7,8 -l 110


    Options:
      -h, --help            show this help message and exit
      -l LEBEDEV, --lebedev=LEBEDEV
                            The number of grid points for the spherical averaging.
                            [default=350]. Select from: 6, 14, 26, 38, 50, 74, 86,
                            110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770,
                            974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470,
                            3890, 4334, 4802, 5294, 5810
      --rlow=RLOW           The smallest radius for the density profile (in
                            angstroms). [default=2e-05]
      --rhigh=RHIGH         The largest radius for the density profile (in
                            angstroms). [default=20.0]
      --num-steps=NUM_STEPS
                            The number of steps in density profile. [default=100]
      --max-ion=MAX_ION     The maximum ionization to consider. [default=2]
      --qc                  Specify the qc convergence scheme in gaussian input.
                            [default=False]




Atomic charges
--------------

Atomic dipoles
--------------

Testing charges and dipoles on the ESP grid
-------------------------------------------

Atomic multipole expansions
---------------------------

Gross and net charges
---------------------

Bond orders, valences and free valences
---------------------------------------

Spin charges
------------

Overlap matrices
----------------

