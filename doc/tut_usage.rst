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

This section is supposed to be read like a tutorial, i.e. in order.


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


Partitioning tools
------------------

All other scripts besides ``hi-atomdb.py`` have the same usage description::

    toon@poony# hi-some-script.py --help
    Usage: hi-some-script.py [options] gaussian.fchk scheme [scheme parameters]
    ...

In the first subsection the usage will be discussed extensively for the script
``hi-charges.py``, but this discussion also applies to all subsequent scripts.
Unless otherwise notices, the example wavefunction is obtained with a HF/STO-3G
computation on hydrogen fluoride in Gaussian03. The example with the spin
charges uses a wavefunction of the OH radical at the same level of theory.

Atomic charges
^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^

Atomic dipoles are compute with the program ``hi-dipoles.py``. Like most HiPart
scripts it takes the same arguments and options as the script ``hi-charges.py``,
which are discussed in the previous section. The online help starts as follows::

    toon@poony# hi-dipoles.py --help
    Usage: hi-dipoles.py [options] gaussian.fchk scheme [scheme parameters]

    hi-dipoles.py computes atomic charges and dipoles.

    These atomic charges and dipoles are the monopole and dipole terms in the
    multipole expansion of each atomic contribution to the density. The atomic
    densities are obtained from the 'scheme' specified at the command line.
    ...

The screen output is also very similar. Depending on the previously executed
scripts, e.g. ``hi-charges.py``, some intermediate results can be loaded from
the work directory and do not have to be computed again.

The dipoles in the Hirshfeld-I scheme can for example be computed as follows::

    toon@poony# hi-dipoles.py gaussian.fchk hirshi atoms/densities.txt
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
    BEGIN Atomic dipoles
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
        BEGIN Iterative Hirshfeld
          Iteration 000    max change = 1.37747e-01    total charge = -3.17836e-05
          Iteration 001    max change = 3.95107e-02    total charge = -3.50733e-05
          Iteration 002    max change = 1.19662e-02    total charge = -3.55995e-05
          Iteration 003    max change = 3.68860e-03    total charge = -3.57148e-05
          Iteration 004    max change = 1.14334e-03    total charge = -3.57459e-05
          Iteration 005    max change = 3.55013e-04    total charge = -3.57552e-05
          Iteration 006    max change = 1.10292e-04    total charge = -3.57580e-05
          Iteration 007    max change = 3.42705e-05    total charge = -3.57589e-05
        END Iterative Hirshfeld
        Writing results to workdir
      END Defining atomic weight functions (each on their own atomic grid)
      Computing dipoles:
       0% ... 100%
      Written gaussian.hipart/hirshi_dipoles.txt
    END Atomic dipoles
    toon@poony# cat gaussian.hipart/hirshi_dipoles.txt
    number of atoms: 2
      i        Z      Dipole-X        Dipole-Y        Dipole-Z      Dipole-norm
    -------------------------------------------------------------------------------
      1   F    9  -0.000007002399  0.000004713804 -0.069590915861  0.069590916373
      2   H    1   0.000002672409  0.000007186970 -0.033744282116  0.033744282987


Testing charges and dipoles on the ESP grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is well known that some charge population schemes are better in reproducing
the electrostatic potential around a molecule than others. [Bultinck2009]_ Since this
often a desirable property, the script ``hi-esp-test.py`` can be used to verify
the quality of the atomic charges and/or the dipoles obtained with a
partitioning scheme. This script computes dipole moment based on charges,
dipoles, and charges with dipoles and compares it with the dipole vector
reported in the formatted checkpoint file. A second test is based on an ESP
fitting cost function for the charges and the dipoles. With this cost function
the script computes how well the charges, the dipoles, and the charges with
the dipoles reproduce the ESP around the molecule.

The current definition of the ESP cost function is discussed in
[Verstraelen2009]_, and we quote the relevant paragraph below for the details:

    We do not rely on charges that are fitted to reproduce the
    ESP around the molecule because they generally suffer from
    statistical inaccuracies. This does not mean that the ESP
    around the molecule is an irrelevant quantity. For the development
    of the electrostatic term in a FF model, one is, in
    principle, only interested in the reproduction of the ESP generated
    by the full electron density, not only in the gas phase
    but also when the electron density adapts to an electrostatic
    perturbation. Under these conditions one can reproduce the
    correct electrostatic interactions. We evaluated, for each
    single point calculation, the ab initio ESP on a molecular
    grid to benchmark the performance of each parametrization.
    A two-dimensional schematic picture of the grid is given in
    Fig. 4. It is constructed as follows. First, 30 concentric
    spheres are placed around each atom. The minimum sphere
    radius is 1.5 times the radius of the noble gas core of the
    corresponding atom, the maximum radius is 30 times the
    noble gas core radius. The radii of intermediate spheres are
    equidistant on a logarithmic scale. On each sphere, we used
    randomly rotated 50-point Lebedev–Laikov grids. The
    random rotation avoids arbitrary preferred directions. For
    this study, we only retained the grid points where the electron
    density is lower than 10e−5 a.u.

This is figure 4 from the paper:

.. image:: grid.png

Again, the script is executed in the same style as all other scripts. See the
documentation of ``hi-charges.py`` for more details. The example below tests the
charges and dipoles obtained with a regular Hirshfeld partitioning::

    toon@poony# hi-esp-test.py gaussian.fchk hirsh atoms/densities.txt
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
    BEGIN Testing charges and dipoles on ESP grid.
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
          BEGIN Conventional Hirshfeld (with neutral pro-atoms)
          END Conventional Hirshfeld (with neutral pro-atoms)
          Writing results to workdir
        END Defining atomic weight functions (each on their own atomic grid)
        Computing charges:
         0% ... 100%
        Written gaussian.hipart/hirsh_charges.txt
      END Atomic charges
      BEGIN Atomic dipoles
        Computing dipoles:
         0% ... 100%
        Written gaussian.hipart/hirsh_dipoles.txt
      END Atomic dipoles
      BEGIN Computing the ESP cost function
        BEGIN Molecular density on the molecular grid
          BEGIN Molecular grid
            BEGIN Estimating noble gas core radii
              Computing noble radii
            END Estimating noble gas core radii
            Constructing molecular grid:
             0% ........... 33% .......... 66% .......... 100%
          END Molecular grid
        END Molecular density on the molecular grid
        BEGIN Molecular potential on the molecular grid
          This may take a minute. Hang on.
        END Molecular potential on the molecular grid
        Written gaussian.hipart/mol_esp_cost.txt
      END Computing the ESP cost function
      Written gaussian.hipart/hirsh_esp_test.txt
    END Testing charges and dipoles on ESP grid.
    toon@poony # ls gaussian.hipart
    hirsh_charges.txt
    hirsh_dipoles.txt
    hirsh_esp_test.txt
    mol_esp_cost.txt
    work

This script computes the charges and dipoles with the given scheme if they are
not present yet. Then the matrix representation of the cost function is
constructed and stored in the file ``mol_esp_cost.txt``. The results of the
test are written in ``hirsh_esp_test.txt``. The output in this example is::

    Reproduction of the molecular dipole
    -------------------------------------------------------------------------------
                      Dipole-X        Dipole-Y        Dipole-Z       Dipole-norm
    -------------------------------------------------------------------------------
    charges (q)    0.000000000000  0.000000000000 -0.262335684109  0.262335684109
    dipoles (p)    0.000001835135  0.000007828203 -0.211588743075  0.211588743228
    q and p        0.000001835135  0.000007828203 -0.473924427184  0.473924427252
    total density  0.000000000000  0.000000000000 -0.473896291000  0.473896291000
    -------------------------------------------------------------------------------

    Reproduction of the external molecular ESP
    -------------------------------------------------------------
                         RMSD             RMS       CORRELATION
    -------------------------------------------------------------
    charges (q)       3.50958e-03     4.39779e-03       1.00
    dipoles (p)       4.36186e-03     3.60738e-03       0.99
    q and p           9.25998e-04     7.99700e-03       0.99
    total density                     7.87933e-03
    -------------------------------------------------------------

As can be seen in the first section, the charges with the dipoles are able to
reproduce the QM dipole moment from Gaussian up to some numerical error. This
error can be controlled to some extent by tuning the grids. In principle, the
correspondence should be exact.

In the second section of the output the QM ESP on the grid points is compared
with the ESP generated by either the charges, the dipoles or the charges with
dipoles. The first column is the root means square deviation over all grid
points. The second column contains the root means square value of the ESP over
all grid points. The third column contains the correlation coefficient between
the approximate and QM ESP data.

In this example it is clear that charges combined with dipoles give already a
fairly accurate description of the ESP around the molecule.

Atomic multipole expansions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The multipole expansion of each atom, up to the hexadecapole, is computed with
the script ``hi-multipoles.py``. The multipoles in the output are computed using
the following real solid harmonics:

========   ======================================================================================================================================================
   (0,0)   :math:`1`
   (1,0)   :math:`z`
  (1,1+)   :math:`x`
  (1,1-)   :math:`y`
   (2,0)   :math:`z^{2} - \frac{1}{2} x^{2} - \frac{1}{2} y^{2}`
  (2,1+)   :math:`x z \sqrt{3}`
  (2,1-)   :math:`y z \sqrt{3}`
  (2,2+)   :math:`\frac{1}{2} \sqrt{3} x^{2} - \frac{1}{2} \sqrt{3} y^{2}`
  (2,2-)   :math:`x y \sqrt{3}`
   (3,0)   :math:`- \frac{3}{2} z x^{2} - \frac{3}{2} z y^{2} + z^{3}`
  (3,1+)   :math:`x \sqrt{6} z^{2} - \frac{1}{4} x \sqrt{6} y^{2} - \frac{1}{4} \sqrt{6} x^{3}`
  (3,1-)   :math:`y \sqrt{6} z^{2} - \frac{1}{4} y \sqrt{6} x^{2} - \frac{1}{4} \sqrt{6} y^{3}`
  (3,2+)   :math:`\frac{1}{2} z \sqrt{15} x^{2} - \frac{1}{2} z \sqrt{15} y^{2}`
  (3,2-)   :math:`x y z \sqrt{15}`
  (3,3+)   :math:`- \frac{3}{4} x \sqrt{10} y^{2} + \frac{1}{4} \sqrt{10} x^{3}`
  (3,3-)   :math:`\frac{3}{4} y \sqrt{10} x^{2} - \frac{1}{4} \sqrt{10} y^{3}`
   (4,0)   :math:`- 3 x^{2} z^{2} - 3 y^{2} z^{2} + \frac{3}{4} x^{2} y^{2} + z^{4} + \frac{3}{8} x^{4} + \frac{3}{8} y^{4}`
  (4,1+)   :math:`- \frac{3}{4} x z \sqrt{10} y^{2} + x \sqrt{10} z^{3} - \frac{3}{4} z \sqrt{10} x^{3}`
  (4,1-)   :math:`- \frac{3}{4} y z \sqrt{10} x^{2} + y \sqrt{10} z^{3} - \frac{3}{4} z \sqrt{10} y^{3}`
  (4,2+)   :math:`- \frac{3}{2} \sqrt{5} y^{2} z^{2} + \frac{3}{2} \sqrt{5} x^{2} z^{2} - \frac{1}{4} \sqrt{5} x^{4} + \frac{1}{4} \sqrt{5} y^{4}`
  (4,2-)   :math:`3 x y \sqrt{5} z^{2} - \frac{1}{2} x \sqrt{5} y^{3} - \frac{1}{2} y \sqrt{5} x^{3}`
  (4,3+)   :math:`- \frac{3}{4} x z \sqrt{70} y^{2} + \frac{1}{4} z \sqrt{70} x^{3}`
  (4,3-)   :math:`\frac{3}{4} y z \sqrt{70} x^{2} - \frac{1}{4} z \sqrt{70} y^{3}`
  (4,4+)   :math:`- \frac{3}{4} \sqrt{35} x^{2} y^{2} + \frac{1}{8} \sqrt{35} x^{4} + \frac{1}{8} \sqrt{35} y^{4}`
  (4,4-)   :math:`\frac{1}{2} y \sqrt{35} x^{3} - \frac{1}{2} x \sqrt{35} y^{3}`
========   ======================================================================================================================================================

This is an example with the Iterative Stockholder Analysis::

    toon@poony# hi-multipoles gaussian.fchk isa
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
    BEGIN Atomic multipoles (up to hexadecapols)
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
        BEGIN Iterative Stockholder Analysis
          Generating initial guess for the pro-atoms
          Iteration 000    max change = 1.83586e-01    total charge = -1.74011e-05
          Iteration 001    max change = 1.41107e-02    total charge = 2.20079e-05
          Iteration 002    max change = 5.77635e-03    total charge = 2.30539e-05
          Iteration 003    max change = 2.63903e-03    total charge = 2.23317e-05
          Iteration 004    max change = 1.25533e-03    total charge = 2.23731e-05
          Iteration 005    max change = 6.12593e-04    total charge = 2.23636e-05
          Iteration 006    max change = 3.04249e-04    total charge = 2.23570e-05
          Iteration 007    max change = 1.52854e-04    total charge = 2.23609e-05
          Iteration 008    max change = 7.72321e-05    total charge = 2.24126e-05
        END Iterative Stockholder Analysis
        Writing results to workdir
      END Defining atomic weight functions (each on their own atomic grid)
      Computing multipoles:
       0% ... 100%
      Written gaussian.hipart/isa_multipoles.txt
    END Atomic multipoles (up to hexadecapols)
    toon@poony# cat gaussian.hipart/isa_multipoles.txt
       Multipoles   |      (0,0)           (1,0)           (1,1+)          (1,1-)          (2,0)           (2,1+)          (2,1-)          (2,2+)          (2,2-)          (3,0)           (3,1+)          (3,1-)          (3,2+)          (3,2-)          (3,3+)          (3,3-)          (4,0)           (4,1+)          (4,1-)          (4,2+)          (4,2-)          (4,3+)          (4,3-)          (4,4+)          (4,4-)
    ----------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      1   F    9    | -0.208473542108 -0.056962384547 -0.000003484828 -0.000003288445  0.226609281756  0.000007104614  0.000009067304 -0.000002056780  0.000003654392  0.095657907380  0.000002499200 -0.000015367877 -0.000000510212 -0.000003563595 -0.000019877664 -0.000003674420  0.021043882212 -0.000040662841  0.000035545551  0.000001380811 -0.000039658470  0.000082557663 -0.000008154178 -0.000000116846  0.000038077725
      2   H    1    |  0.208495954720 -0.019739123492 -0.000012588587 -0.000029008410  0.033088605078 -0.000039695523 -0.000107783990  0.000049301682 -0.000008612389  0.025464488897 -0.000099060049 -0.000173129338  0.000130733479  0.000103548935 -0.000140504106  0.000051467706 -0.018290123225 -0.000092722164 -0.000192502479 -0.000215145579  0.000396793174 -0.000724428107  0.000450717665 -0.000336080516  0.000318860082


Net populations
^^^^^^^^^^^^^^^

Net electron populations are computed with the script ``hi-net.py``. This
example computes the net populations using the becke scheme::

    toon@poony# hi-net.py gaussian.fchk becke
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
    BEGIN Net populations
      BEGIN Atomic grids
        Computing/Loading atomic grids (and distances):
         0% ..... 100%
      END Atomic grids
      BEGIN Molecular density on atomic grids
        Computing/Loading densities:
         0% ... 100%
      END Molecular density on atomic grids
      BEGIN Atomic charges
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
      BEGIN Atomic weights on other atoms' grids.
        Computing off-diagonal atom weights:
         0% ..... 100%
      END Atomic weights on other atoms' grids.
      Integrating over products of stockholder weights:
       0% .... 100%
    END Net populations
    toon@poony# cat gaussian.hipart/becke_net_populations.txt
    number of atoms: 2
          Net       |       1  F            2  H
    ----------------+----------------------------------
      1   F    9    |  9.061921975582  0.138512952959
      2   H    1    |  0.138512952959  0.661118433196

The output is a symmetric matrix with net population charges for each atom pair.
The atomic net populations are put on the diagonal, while the bond net
populations are off-diagonal elements.

Bond orders, valences and free valences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The computation of bond orders is currently only supported for SCF computations,
i.e. Hartree Fock and Density Functional Theory. The following example uses
Hirshfeld partitions. Note that the first section of the screen output must
contain ``Orbitals present: True`` for this script to work. ::

    toon@poony# hi-bond-orders.py gaussian.fchk hirsh atoms/densities.txt
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
    BEGIN Bond orders and atomic valences
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
          BEGIN Conventional Hirshfeld (with neutral pro-atoms)
          END Conventional Hirshfeld (with neutral pro-atoms)
          Writing results to workdir
        END Defining atomic weight functions (each on their own atomic grid)
        Computing charges:
         0% ... 100%
        Written gaussian.hipart/hirsh_charges.txt
      END Atomic charges
      BEGIN Atomic overlap matrix elements
        BEGIN Evaluating orbitals on atomic grids
          Computing/Loading orbitals:
           0% ... 100%
        END Evaluating orbitals on atomic grids
        Computing atomic overlap matrices (alpha):
         0% ... 100%
        Written gaussian.hipart/hirsh_alpha_overlap.txt
      END Atomic overlap matrix elements
      Computing bond orders:
       0% .... 100%
      Written gaussian.hipart/hirsh_bond_orders.txt
      Written gaussian.hipart/hirsh_valences.txt
      Written gaussian.hipart/hirsh_free_valences.txt
    END Bond orders and atomic valences

There are three output files::

    toon@poony# cat gaussian.hipart/hirsh_bond_orders.txt
    number of atoms: 2
       Bond order   |       1  F            2  H
    ----------------+----------------------------------
      1   F    9    |  0.000000000000  1.173854539694
      2   H    1    |  1.173854539694  0.000000000000
    toon@poony# cat gaussian.hipart/hirsh_valences.txt
    number of atoms: 2
      i        Z    Valences
    --------------------------------
      1   F    9    1.173773502217
      2   H    1    1.173862471557
    --------------------------------
    toon@poony# cat gaussian.hipart/hirsh_free_valences.txt
    number of atoms: 2
      i        Z  Free valences
    --------------------------------
      1   F    9   -0.000081037477
      2   H    1    0.000007931863
    --------------------------------

The bond orders are written as off-diagonal elements in the first output file.
The diagonal elements are always zero. The atomic valence and free valence are
written in the two following output files, respectively.

Spin charges
^^^^^^^^^^^^

The computation of spin charges only makes sense in the case of open shell
computations. The example below is an analysis of the OH radial with a minimal
basis set, using the Iterative Hirshfeld scheme. Note that the screen output
must contain the line ``Spin density present: True`` in the beginning. If not,
the reported spin charges are always zero. ::

    toon@poony# hi-spin-charges.py gaussian.fchk hirshi densities.txt
    BEGIN Loading Electronic structure
      Data read from: gaussian.fchk ()
      Restricted: True
      Orbitals present: True
      Spin density present: True
      Number of alpha electrons: 5
      Number of beta electrons: 4
      Number of electrons: 9
      Total charge: 0
      Number of atoms: 2
      Chemical formula: OH
    END Loading Electronic structure
    BEGIN Atomic spin charges
      BEGIN Atomic grids
        Computing/Loading atomic grids (and distances):
         0% ..... 100%
      END Atomic grids
      BEGIN Molecular spin density on atomic grids
        Computing/Loading spin densities:
         0% ... 100%
      END Molecular spin density on atomic grids
      BEGIN Defining atomic weight functions (each on their own atomic grid)
        Trying to load weight functions
        Could not load all weight functions from workdir. Computing them...
        BEGIN Iterative Hirshfeld
          BEGIN Molecular density on atomic grids
            Computing/Loading densities:
             0% ... 100%
          END Molecular density on atomic grids
          Iteration 000    max change = 1.12647e-01    total charge = -4.11725e-05
          Iteration 001    max change = 4.09568e-02    total charge = -4.22415e-05
          Iteration 002    max change = 1.55863e-02    total charge = -4.28892e-05
          Iteration 003    max change = 6.04137e-03    total charge = -4.31639e-05
          Iteration 004    max change = 2.35883e-03    total charge = -4.32746e-05
          Iteration 005    max change = 9.23648e-04    total charge = -4.33186e-05
          Iteration 006    max change = 3.62082e-04    total charge = -4.33360e-05
          Iteration 007    max change = 1.42004e-04    total charge = -4.33428e-05
          Iteration 008    max change = 5.57016e-05    total charge = -4.33455e-05
        END Iterative Hirshfeld
        Writing results to workdir
      END Defining atomic weight functions (each on their own atomic grid)
      Computing spin charges:
       0% ... 100%
      Written gaussian.hipart/hirshi_spin_charges.txt
    END Atomic spin charges
    toon@poony# cat gaussian.hipart/hirshi_spin_charges.txt
    number of atoms: 2
      i        Z  Spin charge
    --------------------------------
      1   O    8    0.971635836892
      2   H    1    0.028378790802
    --------------------------------

Obviously, the spin is located on the oxygen atom.

Overlap matrices
^^^^^^^^^^^^^^^^

The overlap matrices are used, amongst other things, for the computation of the
bond orders. If you are just interested in these matrices, use the command
``hi-overlap.py`` to compute them. The example below demonstrates the usage in
with the Iterative Stockholder Analysis::

    toon@poony# hi-overlap.py gaussian.fchk isa
    BEGIN Loading Electronic structure
      Data read from: gaussian.fchk ()
      Restricted: True
      Orbitals present: True
      Spin density present: True
      Number of alpha electrons: 5
      Number of beta electrons: 4
      Number of electrons: 9
      Total charge: 0
      Number of atoms: 2
      Chemical formula: OH
    END Loading Electronic structure
    BEGIN Atomic overlap matrix elements
      BEGIN Atomic grids
        Computing/Loading atomic grids (and distances):
         0% ..... 100%
      END Atomic grids
      BEGIN Evaluating orbitals on atomic grids
        Computing/Loading orbitals:
         0% ... 100%
      END Evaluating orbitals on atomic grids
      BEGIN Defining atomic weight functions (each on their own atomic grid)
        Trying to load weight functions
        Could not load all weight functions from workdir. Computing them...
        BEGIN Iterative Stockholder Analysis
          BEGIN Molecular density on atomic grids
            Computing/Loading densities:
             0% ... 100%
          END Molecular density on atomic grids
          Generating initial guess for the pro-atoms
          Iteration 000    max change = 1.49110e-01    total charge = -4.94398e-05
          Iteration 001    max change = 2.49598e-02    total charge = 5.23430e-06
          Iteration 002    max change = 1.06242e-02    total charge = 6.47169e-06
          Iteration 003    max change = 5.05120e-03    total charge = 5.43768e-06
          Iteration 004    max change = 2.52294e-03    total charge = 5.06534e-06
          Iteration 005    max change = 1.29954e-03    total charge = 4.85936e-06
          Iteration 006    max change = 6.84229e-04    total charge = 4.74549e-06
          Iteration 007    max change = 3.66279e-04    total charge = 4.68984e-06
          Iteration 008    max change = 1.98584e-04    total charge = 4.66978e-06
          Iteration 009    max change = 1.08700e-04    total charge = 4.67515e-06
          Iteration 010    max change = 5.99109e-05    total charge = 4.67318e-06
        END Iterative Stockholder Analysis
        Writing results to workdir
      END Defining atomic weight functions (each on their own atomic grid)
      Computing atomic overlap matrices (alpha):
       0% ... 100%
      Written gaussian.hipart/isa_alpha_overlap.txt
    END Atomic overlap matrix elements
    toon@poony# cat gaussian.hipart/isa_alpha_overlap.txt
    number of orbitals: 6
    number of atoms:  2
    Atom 0: O
     9.9908093969e-01 -8.7530947103e-04 -2.8253400529e-04  6.2466047246e-08 -5.5004377471e-08  2.4741879654e-03
    -8.7530947103e-04  8.7281007511e-01  1.5893418788e-01  2.1649948815e-06 -2.5906448968e-06 -1.7597650185e-01
    -2.8253400529e-04  1.5893418788e-01  7.6517021096e-01 -8.2719361716e-07  1.9583115417e-06  3.0263452744e-01
     6.2466047246e-08  2.1649948815e-06 -8.2719361716e-07  9.7360289392e-01 -1.5657416364e-07 -4.2141605808e-06
    -5.5004377471e-08 -2.5906448968e-06  1.9583115417e-06 -1.5657416364e-07  9.7360087488e-01  2.1485005028e-06
     2.4741879654e-03 -1.7597650185e-01  3.0263452744e-01 -4.2141605808e-06  2.1485005028e-06  5.1125313840e-01
    Atom 1: H
     8.5158691055e-04  9.4679375098e-04  2.9658087818e-04  4.4584601237e-05  1.0693282217e-05 -2.5227778435e-03
     9.4679375098e-04  1.2717571915e-01 -1.5894043847e-01  6.9492351821e-06 -1.5377486484e-06  1.7598211125e-01
     2.9658087818e-04 -1.5894043847e-01  2.3485394004e-01 -2.2053512215e-06  4.6433653859e-06 -3.0261315440e-01
     4.4584601237e-05  6.9492351821e-06 -2.2053512215e-06  2.6442118753e-02  2.7773242269e-05 -9.2883957586e-06
     1.0693282217e-05 -1.5377486484e-06  4.6433653859e-06  2.7773242269e-05  2.6419482875e-02  3.1344263164e-06
    -2.5227778435e-03  1.7598211125e-01 -3.0261315440e-01 -9.2883957586e-06  3.1344263164e-06  4.8880557083e-01

When this script is executed with a restricted wavefunction, only the matrices
for the alpha orbitals are computed and printed.
