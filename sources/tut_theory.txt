.. _theory:

Theoretical background
======================

The summary below is just a description of the methods implemented in HiPart. It
is highly recommended to read the cited papers for a more extensive discussion.

The first section discusses four different schemes to define fuzzy atoms:
:ref:`becke`, :ref:`hirshfeld`, :ref:`hirshfeld-i`, :ref:`isa`. All subsequent
sections give an overview of the quantities that can be derived with HiPart once
fuzzy atoms (i.e. atomic weight functions) are defined.



Fuzzy atom partitioning
-----------------------

One way to partition the molecular density,
:math:`\rho_{\text{mol}}(\mathbf{r})`, into atomic contributions is to define a
weight function, :math:`w_A(\mathbf{r})`, for every atom :math:`A` in the
molecule. The function value is in the range :math:`[0,1]`. The atomic denisty
is then defined as

.. math:: \rho_A(\mathbf{r}) = w_A(\mathbf{r})\rho_{\text{mol}}(\mathbf{r})

The weight functions must satisfy the condition

.. math:: \sum_{A=1}^N w_A(\mathbf{r}) = 1, ~~ \forall \mathbf{r} \in \mathbb{R}^3,
    :label: completeness

where N is the number of atoms, to guarantee that the sum of the atomic
densities gives back the molecular density. There are two classes of weight
functions, one where the weight functions are binary functions with only 0 or 1
as possible function values, and the fuzzy weight functions where all
intermediate values are also allowed. One can think of binary weight functions
as a means to divide the entire molecular volume into mutually exclusive atomic
regions. The resulting atomic densities are completely non-overlapping. The
fuzzy weights lead to atoms that have overlapping densities. All partitioning
methods in HiPart are fuzzy-atom partitioning methods.

It is common practice to introduce for each atom a so-called pro-atomic
function, :math:`\rho_A^{\text{pro}}(\mathbf{r})`, as an auxiliary tool to define the
actual weight functions. The weight functions are derived from the pro-atomic
functions as follows:

.. math:: w_A(\mathbf{r}) = \frac{\rho_A^{\text{pro}}(\mathbf{r})}{\sum_{B=1}^N \rho_B^{\text{pro}}(\mathbf{r})}

The denominator in this expression is called the pro-molecular density. This
definition of the weight function always satisfies condition :eq:`completeness`
for a broad class of pro-atomic functions.

.. _becke:

Becke
^^^^^

Becke [Becke1988]_ proposed a partitioning scheme that was in the first place
meant as an auxiliary tool to divide an integral over the entire molecular
volume into a sum of atomic integrals. Each atomic integral is evaluated on a
spherical grid using numerical techniques. This is far more convenient than
constructing a global molecular integration grid.

The Becke weights are designed to be mathematically as simple as possible, only
using simple polynomials of distances between atoms and grid points. Becke
introduces so called atomic cell functions, :math:`P_A(\mathbf{r})`, which play
exactly the same role as the pro-atomic function. The weights are derived from
the cell functions as follows:

.. math:: w_A(\mathbf{r}) = \frac{P_A(\mathbf{r})}{\sum_{B=1}^N P_B(\mathbf{r})}

Each cell function is the product of a series of switching functions
:math:`s_{AB}(\mathbf{r})`:

.. math:: P_A(\mathbf{r}) = \prod_{B\neq A} s_{AB}(\mathbf{r})

The switching function goes smoothly from 1 to 0 as one moves from atom
:math:`A` to atom :math:`B`. The gradient of the switching function becomes zero
in the vicinity of the nuclei. Now it is only a matter of constructing a simple
switching function to complete the partitioning method.

Becke proposed the following approach. He starts from one of the elliptical
coordinates:

.. math:: \mu_{AB}(\mathbf{r}) = \frac{|\mathbf{r}_A - \mathbf{r}| - |\mathbf{r}_B - \mathbf{r}|}{|\mathbf{r}_B - \mathbf{r}_A|},

where :math:`\mathbf{r}_A` is the position of atom :math:`A` and
:math:`\mathbf{r}_B` is the position of atom :math:`B`. This coordinate is -1 at
the position of atom :math:`A` and +1 at the position of atom :math:`B`. Then he
introduces the functions :math:`f_k` as follows:

.. math::
    f_1(x) = \frac{x}{2}(3-x^2)

    f_k(x) = f_1(f_{k-1}(x))

The nice property of these functions :math:`f_k` is that :math:`f_k(\mu_{AB})`
is -1 and +1 at the respective nuclei and that the gradient of this function
becomes zero in the vicinity of the nuclei. A simple transformation of the
function :math:`f_3` is used to define the switching function.

.. math:: s_{AB}(\mathbf{r}) = \frac{1}{2}(1-f_3(\mu_{AB}(\mathbf{r})))

The choice of iteration order, :math:`k`, is somewhat arbitrary, but Becke
experienced that 3 was a good trade-off between the sphericity of the atomic
densities (to limit the density of the integration grids) and the locality of
the atomic densities (to limit the extent of the integration grids).

The above definition of the switching functions (and hence weight functions) is
suitable for homonuclear systems. However, for heteronuclear functions it is
desirable to transform the elliptical coordinate, :math:`\mu_{AB}`, such that it
crosses zero around the point where the density in the bond region has a
saddle point. Becke proposes the following transformation:

.. math::
    \nu_{AB} = \mu_{AB} + a_{AB}(1 - \mu_{AB}^2)
    :label: transform_hetero

    s_{AB,\text{het}}(\mu_{AB}) = s_{AB}(\nu_{AB})

The parameter :math:`a_{AB}` controls the position between atoms :math:`A` and
:math:`B` where :math:`\nu_{AB}` goes through zero, and can be used to tune
the size of the basins defined by the weight functions. Based on the covalent
bond radii, :math:`R_A` and :math:`R_B`, Becke defines

.. math::
    u_{AB} = \frac{R_A-R_B}{R_A+R_B}

    a_{AB} = \frac{u_{AB}}{u_{AB}^2-1}

This choice assigns proportionally larger basins to larger atoms in the
molecule, which further improves the convergence of the numerical integrations
over the atomic grids. Note that the absolute value of :math:`a_{AB}` must be
smaller than :math:`\frac{1}{2}` to guarantee that the transform in equation
:eq:`transform_hetero` is monotonous.

In HiPart the parameter :math:`a_{AB}` is constrained to have an absolute value
smaller than 0.45 to suppress pristine behavior. The covalent radii for HiPart
are taken from [Cordero2008]_.

.. _hirshfeld:

Hirshfeld
^^^^^^^^^

Hirshfeld [Hirshfeld1977]_ proposed a partitioning scheme where pro-atomic
densities are derived from computations on neutral atoms by simply averaging the
atomic density over the angular degrees of freedom,

.. math:: \rho_A^{\text{pro}}(|\mathbf{r} - \mathbf{r}_A|) = \rho_A^{\text{pro}}(r) = \int d\Omega \rho_{A,N=Z}^{\text{atom}}(r,\Omega),

where :math:`\Omega` represents the angular degrees of freedom. Prior to the
application of this partitioning scheme one must setup a database of spherically
averaged atomic densities for all elements that are present in the molecule of
interest. For the sake of consistency, this needs to be carried out with the same
level of theory (and basis set) that is used for the molecular computation.

.. _hirshfeld-i:

Iterative Hirshfeld
^^^^^^^^^^^^^^^^^^^

The choice of neutral pro-atoms in the standard Hirshfeld scheme is somewhat
arbitrary. The Iterative Hirshfeld scheme [Bultinck2007]_ is an extension of the
original method, where one seeks for pro-atomic densities that have the same
number of electrons as the atomic partitions in the molecule.

Bultinck et al. introduce a pro-atomic function with an additional parameter,
:math:`N_A`, ie.e the fractional number of electrons in the pro-atomic density.
For integer values of this parameter, the pro-atomic density is just the
spherical average of the corresponding atom in vacuum:

.. math:: \rho_A^{\text{pro}}(r;N_A) = \int d\Omega \rho_{A,N=N_A}^{\text{atom}}(r,\Omega).

For non-integer values of the parameter :math:`N_A`, the pro-atomic density is a
linear interpolation between the two `neighboring` integer-charged atoms:

.. math:: \rho_A^{\text{pro}}(r;N_A) = (\mathrm{ceil}(N_A)-N_A)\rho_A^{\text{pro}}(r;\mathrm{floor}(N_A)) +
                            (N_A-\mathrm{floor}(N_A))\rho_A^{\text{pro}}(r;\mathrm{ceil}(N_A))

The values :math:`N_A` are obtained in an iterative procedure. Initially, they
are all set to zero, and one computes the populations just like in the original
Hirshfeld scheme. In the subsequent iterations the parameters :math:`N_A` are
set to the populations from the previous iteration and one uses these pro-atoms
to compute the population for the next iteration. This is repeated until the
atomic populations converge, i.e. when the maximum absolute value of the
difference in atomic populations between two iterations drops below a predefined
threshold.

Before one can use the Iterative Hirshfeld methods, one must first construct
a database of pro-atomic densities for all the elements in the molecule under
scrutiny. For each element one must compute different charge states.

This scheme is also referred to as `Hirshfeld-I`.


.. _isa:

Iterative Stockholder Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ISA scheme is another extension to the original Hirshfeld method where one
tries to construct spherically symatric pro-atoms that are as similar as
possible to the atomic partitions in the molecule. [Lillestolen2008]_

The initial pro-atoms are constructed by taking the minimal molecular electron
density as a function of the distance from the nucleus. For numerical reasons
this minimal value constrained to be non-zero:

.. math:: \rho_A^{\text{pro},(0)}(|\mathbf{r} - \mathbf{r}_A|) = \max(\epsilon, \min_{\Omega_A} \rho_{\text{mol}}(|\mathbf{r} - \mathbf{r}_A|,\Omega_A))

where :math:`epsilon` is a small positive number and :math:`\Omega_A` are the
angular degrees of freedom of the spherical coordinate system centered at atom
:math:`A`. In each ISA iteration :math:`k`, the new pro-atoms are taken to be
the spherical average of the atomic densities from the previous iteration.

.. math:: \rho_A^{\text{pro},(k+1)}(|\mathbf{r} - \mathbf{r}_A|) = \int d \Omega_A w_A^{(k)}(|\mathbf{r} - \mathbf{r}_A|) \rho_{\text{mol}}(|\mathbf{r} - \mathbf{r}_A|,\Omega_A)

This is again repeated until the atomic populations converge. Note that this
scheme does not depend on a database of atomic densities.


Atomic properties derived from the density
------------------------------------------

In this section we discuss the properties derived from the atomic electron
densities:

.. math:: \rho_A(\mathbf{r}) = w_A(\mathbf{r})\rho_{\text{mol}}(\mathbf{r}),

where :math:`w_A(\mathbf{r})` is the atomic weight function of atom :math:`A`
obtained with some partitioning scheme and :math:`\rho_{\text{mol}}(\mathbf{r})`
is the molecular electron density.

It may be interesting to see how the molecular density is derived from the
density matrix obtained with a quantum chemical ground state computation. First
the basis functions, :math:`f_m` are evaluated in the point :math:`\mathbf{r}`,
where :math:`m` runs from 1 to the number of basis functions. In matrix notation
the molecular density is then computed as follows:

.. math:: \rho_{\text{mol}}(\mathbf{r}) = (f(\mathbf{r}))^T D f(\mathbf{r}),

where :math:`D` is the density matrix.

Charges, Dipoles & Multipoles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The multipole expansion of an atomic density with respect to its nucleus is
defined as follows:

.. math:: Q_l^m = \int (-\rho_A(\mathbf{r})) R_l^m(\mathbf{r}-\mathbf{r}_A) d\mathbf{r}
    :label: multipole

where :math:`R_l^m` is the regular solid harmonic. This multipole expansion
can be used to compute the external electrostatic potential generated by the
atomic density. By external we mean `outside` the atom where the atomic electron
density is negligible. This potential is computed as follows:

.. math:: V^{\text{ext}}_A(\mathbf{r}) = \sum_{l=0}^{\infty} \sum_{m=-l}^{l} (-1)^m I_l^{-m}(\mathbf{r}-\mathbf{r}_A) Q_l^m

where :math:`I_l^m` is the irregular solid harmonic. The regular and
irregular solid harmonics are related to the spherical harmonics, :math:`Y_l^m`
as follows:

.. math::
    R_l^m(\mathbf{r}) = \sqrt{\frac{4\pi}{2l+1}} r^l Y_l^m(\Omega)

    I_l^m(\mathbf{r}) = \sqrt{\frac{4\pi}{2l+1}} r^{-l-1} Y_l^m(\Omega)

In HiPart we use the real valued variants of these functions and we replace the
integrals :eq:`multipole` by their real counterparts:

.. math:: \tilde{Q}_l^m = \int (-\rho_A(\mathbf{r})) \tilde{R}_l^m(\mathbf{r}-\mathbf{r}_A) d\mathbf{r},
    :label: multipole-real

with

.. math:: \tilde{R}_l^m = \left\lbrace\begin{array}{ll}
    R_l^0 & \text{if  } m=0\\
    {1\over\sqrt2}\left(R_l^m+(-1)^m \, R_l^{-m}\right)  & \text{if  } m>0 \\
    {1\over i\sqrt2}\left(R_l^{-m}-(-1)^{m}\, R_l^{m}\right) & \mbox{if  } m<0.
    \end{array} \right.

The table below lists all regular solid harmonics implemented in HiPart, i.e. up
to the hexadecapole, :math:`(l=4)`. The formulae are automatically generated,
which causes a somewhat ugly formatting.

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

The first line in this table corresponds to the atomic population, after adding
the nuclear charge, one obtains the `effective` atomic charge. The subsequent
three rows correspond to the components of the atomic dipole, and so on.

Net and overlap populations
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The net and overlap populations are obtained by inserting the weight function
twice in the integral over the molecular density:

.. math:: N^{\text{no}}_{AB} = \int w_A(\mathbf{r}) w_B(\mathbf{r}) \rho_{\text{mol}}(\mathbf{r}) d\mathbf{r}

The diagonal elements are the net populations. One can interpret the net
population of an atom as the amount of electrons that is associated only with
that atom. The overlap populations, i.e. off-diagonal elements, can be
interpreted as the amount of electrons that are shared between two atoms.

These quantities are certainly not convenient as measures for atomic valence and
bond order, although one can expect that there must be some correlation between
the overlap population and the bond order. The main issue is that these numbers
are not even close to the integer values that one would expect for atomic
valences and bond orders. In the case of binary weight functions, the overlap
populations would be zero and the net charges would be the regular atomic
populations.


Atomic properties derived from the spin density
-----------------------------------------------

In this section we discuss the properties derived from the atomic spin
densities:

.. math:: \rho^{\text{spin}}_A(\mathbf{r}) = w_A(\mathbf{r})\rho^{\text{spin}}_{\text{mol}}(\mathbf{r}),

where :math:`w_A(\mathbf{r})` is the atomic weight function of atom :math:`A`
obtained with some partitioning scheme and
:math:`\rho^{\text{spin}}_{\text{mol}}(\mathbf{r})` is the molecular spin
density.

The spin density is derived from the spin density matrix in the same way as the
conventional electron density is derived from the density matrix:

.. math:: \rho^{\text{spin}}_{\text{mol}}(\mathbf{r}) = (f(\mathbf{r}))^T D^{\text{spin}} f(\mathbf{r}),

where :math:`f(\mathbf{r})` is the vector with basis functions evaluated in
point :math:`\mathbf{r}` and :math:`D^{\text{spin}}` is the spin density matrix.
The spin density matrix and the conventional density matrix can be derived from
the alpha spin density matrix, :math:`D^{\alpha}`, and the beta spin density
matrix, :math:`D^{\beta}`, as follows:

.. math::
    D^{\text{spin}} = D^{\alpha} - D^{\beta}

    D = D^{\alpha} + D^{\beta}


Spin charges
^^^^^^^^^^^^

The spin charges are the atomic populations derived from the molecular spin
density.

.. math:: N^{\text{spin}}_{A} = \int w_A(\mathbf{r}) \rho^{\text{spin}}_{\text{mol}}(\mathbf{r}) d\mathbf{r}


Atomic overlap matrices (in the basis of contracted Gaussians)
--------------------------------------------------------------

The atomic overlap matrices do not depend on the density or density matrix, but
only depend on the basis functions used to describe the wavefunction. The
conventional overlap matrix is defined as follows:

.. math:: S_{\mu\nu} = \int f_{\mu}(\mathbf{r}) f_{\nu}(\mathbf{r}) d\mathbf{r}

The square root of the overlap matrix can be used to transform the
non-orthogonal basis of contracted Gaussians into an orthonormal basis, and is
in general a tool to work with non-orthogonal basis sets. One defines the atomic
overlap matrix by inserting an atomic weight function into the integral:

.. math:: S^{A}_{\mu\nu} = \int f_{\mu}(\mathbf{r}) w_A(\mathbf{r}) f_{\nu}(\mathbf{r}) d\mathbf{r}

One can define the atomic population as the trace of the product of the density
matrix and the corresponding atomic overlap matrix:

.. math:: N_A = \mathrm{Tr} (D S^A)

Similarly, one can write down the spin charges as function of the spin density
matrix and the overlap matrix:

.. math:: N^{\text{spin}}_A = \mathrm{Tr} (D^{\text{spin}} S^A)


Atomic properties derived from the density matrix
-------------------------------------------------

All properties discussed in the previous sections can be written as an integral
over the electron (spin) density multiplied by a (weight) function.
Alternatively one can also go back to the matrix :math:`(D S^A)` and manipulate
this object before taking the trace. This allows use to derive new quantities
that are not simple written as a function of the electron density.


Bond orders and atomic valences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One defines the bond order, valence and free valence, in the fuzzy atom
framework as follows:

.. math::
    \mathcal{B}_{AB} = 2\mathrm{Tr}\left[(D^{\alpha} S^A)(D^{\alpha} S^B)^T + (D^{\beta} S^A)(D^{\beta} S^B)^T \right]

    \mathcal{V}_A = 2 N_A - \mathrm{Tr}\left[ (D S^A)(D S^A)^T \right]

    \mathcal{F}_A = \mathcal{V}_A - \sum_{B \neq A} \mathcal{B}_{AB}

Mayer has written a `personal account` [Mayer2007]_ about bond orders and
valence indices. It is a good introduction for those who are new to these
concepts.


Atomic overlap matrices (in the basis of the orbitals)
------------------------------------------------------

The overlap matrix can also be computed in the basis of the orbitals. Because
the orbitals are an orthogonal basis, the conventional overlap matrix

.. math:: S^{\text{orb}}_{\mu\nu} = \int \psi_{\mu}(\mathbf{r}) \psi_{\nu}(\mathbf{r}) d\mathbf{r},

where :math:`\psi_{\mu}` and :math:`\psi_{\nu}` are the orbitals, is simply the
identity matrix. One can now define atomic contributions to this matrix by
inserting the atomic weight functions into the integral:

.. math:: S^{\text{orb},A}_{\mu\nu} = \int \psi_{\mu}(\mathbf{r}) w_A(\mathbf{r}) \psi_{\nu}(\mathbf{r}) d\mathbf{r}

This is in principle the same quantity is introduced earlier, but just in a
different basis.

.. _esp:

Electrostatic Potential (ESP) Fitting
-------------------------------------

In this section we will refer to electrostatic potential generated by the
electron density and the nuclei as the full ESP, or :math:`V_{\text{full}}`.

ESP fitting is a procedure where the amplitudes of point monopoles at the
positions of the nuclei are fitted to reproduce the full ESP `outside` the
molecule. Advanced schemes also include dipoles and optionally higher
multipoles, and consider also other critical points than just the positions of
the nuclei.

The goodness of the reproduction of the potential outside the molecule is
typically measured by a fitness (or cost) function. In the most simplistic
approach, this fitness function is simply a sum of weighted squared errors
between the full-blown potential, :math:`V_{\text{full}}`, and the potential
generated by the point charges, :math:`V_{\text{model}}`:

.. math::
    X = \sum_{p=1}^P w_p \left( V_{\text{full}}(\mathbf{r}_p) - V_{\text{model}}(\mathbf{r}_p) \right)^2

    V_{\text{model}}(\mathbf{r}_p) = \sum_{A=1}^N \frac{q_A}{|\mathbf{r}_p - \mathbf{r}_A|}

Minimization of this fitness function with respect to the unknowns, the charges
:math:`q_A`, yields the ESP-optimal charges.

Most ESP fitting methods differ in the way the grid points are constructed. No
matter how one selects the grid points, the cost function :math:`X` is always
ill defined. Several attempts have been made to turn this fitness function into
a well-behaved one, of which the RESP method [Bayly1993]_ is the most
wide-spread. In HiPart, this rank-deficiency issue of the fitness function is less
problematic because HiPart only uses such cost functions to measure how well
charges (and dipoles) derived from a partitioning scheme are able to reproduce
the ESP on a set of grid points around the molecule. HiPart does not compute
ESP-fitted charges.

The selection of grid points for the cost function used in HiPart is discussed
in [Verstraelen2009]_. All weights are set to 1. The relevant paragraph for the
paper is quoted below for the details:

    *We do not rely on charges that are fitted to reproduce the ESP around the
    molecule because they generally suffer from statistical inaccuracies. This
    does not mean that the ESP around the molecule is an irrelevant quantity.
    For the development of the electrostatic term in a FF model, one is, in
    principle, only interested in the reproduction of the ESP generated by the
    full electron density, not only in the gas phase but also when the electron
    density adapts to an electrostatic perturbation. Under these conditions one
    can reproduce the correct electrostatic interactions. We evaluated, for each
    single point calculation, the ab initio ESP on a molecular grid to benchmark
    the performance of each parametrization. A two-dimensional schematic picture
    of the grid is given in Fig. 4. It is constructed as follows. First, 30
    concentric spheres are placed around each atom. The minimum sphere radius is
    1.5 times the radius of the noble gas core of the corresponding atom, the
    maximum radius is 30 times the noble gas core radius. The radii of
    intermediate spheres are equidistant on a logarithmic scale. On each sphere,
    we used randomly rotated 50-point Lebedev–Laikov grids. The random rotation
    avoids arbitrary preferred directions. For this study, we only retained the
    grid points where the electron density is lower than 10e−5 a.u.*

This is figure 4 from the paper:

.. image:: grid.png
