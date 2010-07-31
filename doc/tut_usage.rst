Using HiPart
============

For the terminology used in this chapter we refer to the :ref:`theory`.

HiPart is actually not a single program but a collection of Python scripts that
can be used from the command line to perform a population analysis on a Gaussian
computation as a post-processing step. All scripts start with ``hi-``, which
makes it easy to get an overview of all available script. Just enter ``hi-`` on
the command line and hit the ``TAB`` key twice::

    toon@poony# hi-
    hi-atomdb.py        hi-charges.py       hi-esp-test.py      hi-multipoles.py    hi-spin-charges.py
    hi-bond-orders.py   hi-dipoles.py       hi-gross-net.py     hi-overlap.py
    toon#poony# hi-


Preparing the atomic database
-----------------------------

Both the :ref:`hirshfeld` and :ref:`hirshfeld-i` rely on a database of
pro-atoms. These are spherically averaged atomic densities obtained from
single-atom computations. The program ``hi-atomdb.py`` generates such a database
using the Gaussian03 program. Make sure the g03 program is in the path before
executing ``hi-atomdb.py``.

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
      --qc                  Specify the qc convergence scheme in Gaussian input.
                            [default=False]


The default settings are good enough for most levels of theory and basis sets.
``hi-atomdb.py`` prints some progress information on screen. A typical output
looks like this::

    toon@poony# hi-atomdb.py HF/STO-3G 1,6-10
    Computing atomic database with program Gaussian (g03,qc=False)
    Creating input files:
     0% ......... 33% ......... 66% ......... 100%
    Atomic computations:
     0% ......... 14% ......... 29% ......... 44% ......... 59% ......... 74% ......... 89% ...... 100%
    Density profiles:
     0% ......... 33% ......... 66% ......... 100%
    Total charge error:   1  H -1    1.78596e-09
    Total charge error:   1  H +0    1.68323e-09
    Total charge error:   6  C -2    4.43124e-09
    Total charge error:   6  C -1    4.76319e-09
    Total charge error:   6  C +0    5.58598e-09
    Total charge error:   6  C +1    6.18370e-09
    Total charge error:   6  C +2    6.41107e-09
    Total charge error:   7  N -2    7.81641e-09
    Total charge error:   7  N -1    8.16317e-09
    Total charge error:   7  N +0    8.53874e-09
    Total charge error:   7  N +1    8.58778e-09
    Total charge error:   7  N +2    8.75891e-09
    Total charge error:   8  O -2    -1.60109e-08
    Total charge error:   8  O -1    -1.47804e-08
    Total charge error:   8  O +0    -1.35640e-08
    Total charge error:   8  O +1    -1.25001e-08
    Total charge error:   8  O +2    -1.13051e-08
    Total charge error:   9  F -1    -6.58516e-09
    Total charge error:   9  F +0    -5.80743e-09
    Total charge error:   9  F +1    -4.97888e-09
    Total charge error:   9  F +2    -4.19043e-09
    Total charge error:  10 Ne +0    7.07409e-09
    Total charge error:  10 Ne +1    7.82146e-09
    Total charge error:  10 Ne +2    8.49154e-09

The program consists of three phases: (i) setup of the atomic input files for
Gaussian03, (ii) Gaussian03 computations on every atomic input, and (iii)
derivation of the spherically averaged atomic densities. In the end a check
is performed by integrating the total charge based on the spherically averaged
densities. When the poor grids are used, it will be obvious from the errors in
this last check. In this example the errors are very small because of the
minimal basis set.

The choice of angular grid is not that important an can be chosen very large
because the computations are fast enough anyway. The radial grid settings are
more delicate. The radial grid is always logarithmic, i.e. equidistant on a
logarithmic scale. The same radial grid will be used by all other HiPart
programs that use this database. If for some reason large radial grids are
required later, they have to be defined at this point. In case of Lithium, heavy
atoms or large basis sets, one may want to tune the radial grid.

Once the program is finished, the following files are generated::

    toon@poony# find | sort
    .
    ./001H
    ./001H/neg1
    ./001H/neg1/gs
    ./001H/neg1/mult1
    ./001H/neg1/mult1/gaussian.com
    ./001H/neg1/mult1/gaussian.fchk
    ./001H/neg1/mult1/gaussian.log
    ./001H/neg1/mult1/grid.bin
    ./001H/neg1/mult1/grid_moldens.bin
    ./001H/neg1/mult1/grid_moldens.txt
    ./001H/neg1/mult1/grid.txt

    ...

    ./010Ne/pos2
    ./010Ne/pos2/gs
    ./010Ne/pos2/mult1
    ./010Ne/pos2/mult1/gaussian.com
    ./010Ne/pos2/mult1/gaussian.fchk
    ./010Ne/pos2/mult1/gaussian.log
    ./010Ne/pos2/mult3
    ./010Ne/pos2/mult3/gaussian.com
    ./010Ne/pos2/mult3/gaussian.fchk
    ./010Ne/pos2/mult3/gaussian.log
    ./010Ne/pos2/mult3/grid.bin
    ./010Ne/pos2/mult3/grid_moldens.bin
    ./010Ne/pos2/mult3/grid_moldens.txt
    ./010Ne/pos2/mult3/grid.txt
    ./chieta_au.txt
    ./chieta_ev.txt
    ./densities.txt
    ./energies.txt

For every atom-charge combination, all reasonable spin multiplicities are
computed and the lowest in energy is selected. One can run ``hi-partdb.py`` a
second time with more atoms to extend the database. (Existing computations will
be reused, but make sure the same basis and level of theory are used.)

Only the file ``densities.txt`` will be used later. It has the following format::

    Radii [bohr]               3.7794523e-05 4.3454517e-05 4.9962135e-05 ...
    Densities   1  H -1 [a.u.] 7.8938827e-01 7.8938827e-01 7.8938827e-01 ...
    Densities   1  H +0 [a.u.] 3.9469414e-01 3.9469414e-01 3.9469413e-01 ...
    Densities   6  C -2 [a.u.] 7.9128494e+01 7.9128491e+01 7.9128488e+01 ...
    Densities   6  C -1 [a.u.] 7.9128494e+01 7.9128491e+01 7.9128488e+01 ...
    Densities   6  C +0 [a.u.] 7.9128494e+01 7.9128491e+01 7.9128488e+01 ...
    Densities   6  C +1 [a.u.] 7.9128494e+01 7.9128491e+01 7.9128488e+01 ...
    Densities   6  C +2 [a.u.] 7.9128494e+01 7.9128491e+01 7.9128488e+01 ...
    ...

The first row consists of the radial grid points. All subsequent lines are the
averaged densities of the atom-charge states at the corresponding distances from
the nucleus.


Atomic charges
--------------

Effective atomic charges are computed with ``hi-charges.py``. The online help is
as follows::

    toon@poony# hi-charges.py --help
    Usage: hi-charges.py [options] gaussian.fchk scheme [scheme parameters]

    hi-charges.py computes effective atomic charges.

    The effective atomic charges are the monopole terms in the multipole expansion
    of each atomic contribution to the density plus the monopole of the nucleus. The
    atomic densities are obtained from the 'scheme' specified at the command line.

    The file gaussian.fchk is a formatted checkpoint file from a Gaussian
    computation. To obtain this file, add the following line on top of a Gaussian
    com-file (before running the job)

    %chk=gaussian.chk

    After the Gaussian computation transform this binary checkpoint file into
    a text file with the ``formchk`` program of the Gaussian software suite:

    formchk gaussian.chk gaussian.fchk

    Partitioning schemes:

     * Becke's Smooth Voronoi Partitioning
         scheme = becke
         scheme parameters = [k] [r_low r_high steps]

         The parameter k is optional and defaults to 3. It is the number of
         iterations in the definition of the weight function in Becke's paper.

         Three additional parameters can be provided of the file rs.bin is not yet
         present in the work directory. The first two, r_low and r_high, are the
         first and the last point on the logarithmic radial grid in angstrom. The
         third, steps, is the number of grid points on the radial grid. The default
         is 2.0e-5, 20.0 and 100, respectively.

         Becke, A. D. J. Chem. Phys. 1988,  88, 2547-2553.
         http://dx.doi.org/10.1063/1.454033

     * Hirshfeld Partitioning
         scheme = hirsh
         scheme parameters = densities.txt

         The file densities.txt is generated with the script hi-atomdb.py. It
         contains spherically averaged densities of individual atoms. Make sure all
         the atoms present in the molecule of interest are included in the file
         densities.txt

         Hirshfeld, F. L. Theor. Chim. Acta 1977, 44, 129-138.
         http://dx.doi.org/10.1007/BF00549096

     * Hirshfeld-I Partitioning
         scheme = hirshi
         scheme parameters = densities.txt

         The file densities.txt is generated with the script hi-atomdb.py. It
         contains spherically averaged densities of individual atoms. Make sure all
         the atoms present in the molecule of interest are included in the file
         densities.txt

         Bultinck, P.;  Van Alsenoy, C.;  Ayers, P. W.;  Dorca, R. C. J. Chem. Phys.
         2007, 126, 144111.
         http://dx.doi.org/10.1063/1.2715563

     * Iterative Stockholder Partitioning
         scheme = isa
         scheme parameters = [r_low r_high steps]

         Three additional parameters can be provided of the file rs.bin is not yet
         present in the work directory. The first two, r_low and r_high, are the
         first and the last point on the logarithmic radial grid in angstrom. The
         third, steps, is the number of grid points on the radial grid. The default
         is 2.0e-5, 20.0 and 100, respectively.

         Lillestolen, T. C.;  Wheatley, R. J. Chem. Commun. 2008,  5909-5911.
         http://dx.doi.org/10.1039/b812691g



    Options:
      -h, --help            show this help message and exit
      -l LEBEDEV, --lebedev=LEBEDEV
                            The number of grid points for the atomic grids.
                            [default=110]. Select from: 6, 14, 26, 38, 50, 74, 86,
                            110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770,
                            974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470,
                            3890, 4334, 4802, 5294, 5810
      -c, --clean           If given, the workdir with the binary data is removed
                            when the computation has finished.
      -n, --no-fix-total-charge
                            Do not correct the total charge.
      -t THRESHOLD, --threshold=THRESHOLD
                            When the maximum change in the charges drops below
                            this threshold value, the iteration stops.
                            [default=0.0001]
      --max-iter=MAX_ITER   Maximum number of iterations in self-consistent
                            procedures. [default=500]


Note that there are four different Fuzzy atom schemes that can be used to define
atomic populations. The screen output gives some progress information of the
various parts of the program. This is an example screen output::

    toon@poony# hi-charges.py gaussian.fchk becke
    BEGIN Loading Electronic structure
      Data read from: gaussian.fchk ()
      Restricted: True
      Orbitals present: True
      Spin density present: False
      Number of alpha electrons: 5
      Number of beta electrons: 5
      Number of electrons: 10
      Total charge: 0
      Number of atoms: 2
      Chemical formula: FH
    END Loading Electronic structure
    BEGIN Atomic charges
      BEGIN Atomic grids
        Computing/Loading atomic grids (and distances):
         0% ..... 100%
      END Atomic grids
      BEGIN Molecular density on atomic grids
        Computing/Loading densities:
         0% ... 100%
      END Molecular density on atomic grids
      BEGIN Defining atomic weight functions (each on their own atomic grid)
        Trying to load weight functions
        Could not load all weight functions from workdir. Computing them...
        BEGIN Becke's Smooth Voronoi Partitioning
          Computing/Loading cell functions:
           0% ... 100%
        END Becke's Smooth Voronoi Partitioning
        Writing results to workdir
      END Defining atomic weight functions (each on their own atomic grid)
      Computing charges:
       0% ... 100%
      Written gaussian.hipart/becke_charges.txt
    END Atomic charges


The entire screen output is conceived as a call graph that shows in which part
of Hipart the program is currently active. The order of the routines is
determined by an internal dependency resolver that allows many different
workflows through the program. The first part of the output is a summary of the
electronic structure stored in the file ``gaussian.fchk``. From then on the
actual computation is carried out.

All output is stored in a subdirectory of the current directory whose name is
based on the filename of the formatted checkpoint file. E.g. if the formatted
checkpoint file is ``gaussian.fchk``, then the output directory is
``gaussian.hipart``. In this example the following output files can be found in
``gaussian.hipart``::

    toon@poony# ls gaussian.hipart/
    becke_charges.txt
    work

All output that depends on the choice of the partitioning scheme is prefixed
with the corresponding key, e.g. in this case we have ``becke_charges.txt``. The
work directory contains cached binary intermediate results that will be reused
when another HiPart script (or the same script with different options) is
executed afterwards. It can always be removed, or with the ``--clean`` option it
is automatically removed. In this example the ``work`` directory contains the
following files::

    toon@poony# ls gaussian.hipart/work/
    atom00000_becke_atweights.bin  atom00001_cell00000.bin
    atom00000.bin                  atom00001_cell00001.bin
    atom00000_cell00000.bin        atom00001_moldens.bin
    atom00000_cell00001.bin        becke_charges.bin
    atom00000_moldens.bin          becke_populations.bin
    atom00001_becke_atweights.bin  context
    atom00001.bin                  rs.bin

Certain choices (grids and some other options) affect the content of the files
in the work directory. When different grids are used in a second run, the work
directory is no longer usable and you will get an error message like this::

    toon@poony# hi-charges.py gaussian.fchk becke -l14
    BEGIN Electronic structure summary
      Data read from: gaussian.fchk ()
      Restricted: True
      Orbitals present: True
      Spin density present: False
      Number of alpha electrons: 5
      Number of beta electrons: 5
      Number of electrons: 10
      Total charge: 0
      Number of atoms: 2
      Chemical formula: FH
    END Electronic structure summary
    Traceback (most recent call last):
      File "/home/toon/bin/hi-charges.py", line 5, in <module>
        pkg_resources.run_script('HiPart==0.004', 'hi-charges.py')
      File "/usr/lib/python2.6/dist-packages/pkg_resources.py", line 461, in run_script
        self.require(requires)[0].run_script(script_name, ns)
      File "/usr/lib/python2.6/dist-packages/pkg_resources.py", line 1194, in run_script
        execfile(script_filename, namespace, namespace)
      File "/home/toon/lib/python/HiPart-0.004-py2.6-linux-x86_64.egg/EGG-INFO/scripts/hi-charges.py", line 35, in <module>
        context, cache = parse_command_line(usage)
      File "/home/toon/lib/python/HiPart-0.004-py2.6-linux-x86_64.egg/hipart/opts.py", line 98, in parse_command_line
        cache = CacheClass.new_from_args(context, args[2:])
      File "/home/toon/lib/python/HiPart-0.004-py2.6-linux-x86_64.egg/hipart/cache.py", line 1046, in new_from_args
        return cls(context, k, rs)
      File "/home/toon/lib/python/HiPart-0.004-py2.6-linux-x86_64.egg/hipart/cache.py", line 1051, in __init__
        BaseCache.__init__(self, context, {"becke_k": str(k)})
      File "/home/toon/lib/python/HiPart-0.004-py2.6-linux-x86_64.egg/hipart/cache.py", line 99, in __init__
        self.context.check_tag(extra_tag_attributes)
      File "/home/toon/lib/python/HiPart-0.004-py2.6-linux-x86_64.egg/hipart/context.py", line 85, in check_tag
        raise ContextError("The existing work directory contains incompatible data. Trash it!")
    hipart.context.ContextError: The existing work directory contains incompatible data. Trash it!

Either remove the entire work directory, or stick to the options used in the
first execution of a HiPart script.

The output file ``becke_charges.txt`` has the following contents::

    number of atoms: 2
      i        Z      Charge
    --------------------------------
      1   F    9   -0.200384125318
      2   H    1    0.200384125318
    --------------------------------

It is easily processed with other programs in a follow-up analysis. Note that
the same data are also present in binary format in the file
``work/becke_charges.bin.``

Although the numbers in the output file are printed with 13 decimals, one must
realize that precision is not the same as accuracy. The accuracy of these
numbers depends on the choice of the radial and angular grids. The accuracy is
also inherently limited by the choices made in the Gaussian input file and the
precision of the numbers in the formatted checkpoint file.

Even a second run of the program (after removing the work directory) will result
in slightly different numbers::

    toon@poony# rm -r gaussian.hipart
    toon@poony# hi-charges.py gaussian.fchk becke
    toon@poony# cat gaussian.hipart/becke_charges.txt
    number of atoms: 2
      i        Z      Charge
    --------------------------------
      1   F    9   -0.200392380104
      2   H    1    0.200392380104
    --------------------------------

This is due to the random rotations applied to the angular grids. This practice
slightly has several advantages:

* It improves the accuracy due to compensation of errors.
* It removes directionional preference in the grids and.
* It allows simple estimates of the accuracy by simply rerunning the same
  analysis twice.

For the sake of completeness, these are the commands to compute the charges on
the same molecule with the three other partitioning schemes::

    toon@poony# hi-charges.py gaussian.fchk hirsh atoms/densities.txt
    toon@poony# hi-charges.py gaussian.fchk hirshi atoms/densities.txt
    toon@poony# hi-charges.py gaussian.fchk isa


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
