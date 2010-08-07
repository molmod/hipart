.. _theory:

Theoretical background
======================

The summary below is just a description of the methods implemented in HiPart. It
is highly recommended to read the cited papers for a more extensive discussion.

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


Atomic properties derived from the partitioned density
------------------------------------------------------

In this section we discuss the properties derived from the atomic electron
densities:

.. math:: \rho_A(\mathbf{r}) = w_A(\mathbf{r})\rho_{\text{mol}}(\mathbf{r}),

where :math:`w_A(\mathbf{r})` is the atomic weight function of atom :math:`A`
obtained with some partitioning scheme and :math:`\rho_{\text{mol}}` is the
molecular electron density.

It may be interesting to see how the molecular density is derived from the
density matrix obtained with a quantum chemical ground state computation. First
the basis functions, :math:`B_m` are evaluated in the point :math:`\mathbf{r}`,
where :math:`m` runs from 1 to the number of basis functions. In matrix notation
the molecular density is then computed as follows:

.. math:: \rho_{\text{mol}}(\mathbf{r}) = (B(\mathbf{r}))^T D B(\mathbf{r}),

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


Atomic properties derived from the partitioned spin density
-----------------------------------------------------------

TODO

Spin charges
^^^^^^^^^^^^

TODO

Atomic properties derived from the partitioned density matrix
-------------------------------------------------------------

TODO

Bond orders
^^^^^^^^^^^

TODO

Atomic overlap matrices (in the basis of contracted Gaussians)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Atomic properties derived from the partitioned orbitals
-------------------------------------------------------

Atomic overlap matrices (in the basis of the orbitals)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Electrostatic Potential Fitting
-------------------------------

TODO
