.. _theory:

Theoretical background
======================

The summary below is just a description of the methods implemented in HiPart. It
is highly recommended to read the cited papers for a more extensive discussion.

Fuzzy atom partitioning
-----------------------

One way to partition the molecular density,
:math:`\rho_{\text{mol}}(\mathbf{r})`, into atomic contributions is to define a
weight function, :math:`w_i(\mathbf{r})`, for every atom :math:`i` in the
molecule. The function value is in the range :math:`[0,1]`. The atomic denisty
is then defined as

.. math:: \rho_i(\mathbf{r}) = w_i(\mathbf{r})\rho_{\text{mol}}(\mathbf{r})

The weight functions must satisfy the condition

.. math:: \sum_{i=1}^N w_i(\mathbf{r}) = 1, ~~ \forall \mathbf{r} \in \mathbb{R}^3,
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
function, :math:`\rho_i^{\text{pro}}(\mathbf{r})`, as an auxiliary tool to define the
actual weight functions. The weight functions are derived from the pro-atomic
functions as follows:

.. math:: w_i(\mathbf{r}) = \frac{\rho_i^{\text{pro}}(\mathbf{r})}{\sum_{j=1}^N \rho_j^{\text{pro}}(\mathbf{r})}

The denominator in this expression is called the pro-molecular density. This
definition of the weight function always satisfies condition :eq:`completeness`
for a broad class of pro-atomic functions.

.. _becke:

Becke scheme
^^^^^^^^^^^^

Becke [Becke1988]_ proposed a partitioning scheme that was in the first place
meant as an auxiliary tool to divide an integral over the entire molecular
volume into a sum of atomic integrals. Each atomic integral is evaluated on a
spherical grid using numerical techniques. This is far more convenient than
constructing a global molecular integration grid.

The Becke weights are designed to be mathematically as simple as possible, only
using simple polynomials of distances between atoms and grid points. Becke
introduces so called atomic cell functions, :math:`P_i(\mathbf{r})`, which play
exactly the same role as the pro-atomic function. The weights are derived from
the cell functions as follows:

.. math:: w_i(\mathbf{r}) = \frac{P_i(\mathbf{r})}{\sum_{j=1}^N P_j(\mathbf{r})}

Each cell function is the product of a series of switching functions
:math:`s_{ij}(\mathbf{r})`:

.. math:: P_i(\mathbf{r}) = \prod_{j\neq i} s_{ij}(\mathbf{r})

The switching function goes smoothly from 1 to 0 as one moves from atom
:math:`i` to atom :math:`j`. The gradient of the switching function becomes zero
in the vicinity of the nuclei. Now it is only a matter of constructing a simple
switching function to complete the partitioning method.

Becke proposed the following approach. He starts from one of the elliptical
coordinates:

.. math:: \mu_{ij}(\mathbf{r}) = \frac{|\mathbf{r}_i - \mathbf{r}| - |\mathbf{r}_j - \mathbf{r}|}{|\mathbf{r}_j - \mathbf{r}_i|},

where :math:`\mathbf{r}_i` is the position of atom :math:`i` and
:math:`\mathbf{r}_j` is the position of atom :math:`j`. This coordinate is -1 at
the position of atom :math:`i` and +1 at the position of atom :math:`j`. Then he
introduces the functions :math:`f_k` as follows:

.. math::
    f_1(x) = \frac{x}{2}(3-x^2)

    f_k(x) = f_1(f_{k-1}(x))

The nice property of these functions :math:`f_k` is that :math:`f_k(\mu_{ij})`
is -1 and +1 at the respective nuclei and that the gradient of this function
becomes zero in the vicinity of the nuclei. A simple transformation of the
function :math:`f_3` is used to define the switching function.

.. math:: s_{ij}(\mathbf{r}) = \frac{1}{2}(1-f_3(\mu_{ij}(\mathbf{r})))

The choice of iteration order, :math:`k`, is somewhat arbitrary, but Becke
experienced that 3 was a good trade-off between the sphericity of the atomic
densities (to limit the density of the integration grids) and the locality of
the atomic densities (to limit the extent of the integration grids).

The above definition of the switching functions (and hence weight functions) is
suitable for homonuclear systems. However, for heteronuclear functions it is
desirable to transform the elliptical coordinate, :math:`\mu_{ij}`, such that it
crosses zero around the point where the density in the bond region has a
saddle point. Becke proposes the following transformation:

.. math::
    \nu_{ij} = \mu_{ij} + a_{ij}(1 - \mu_{ij}^2)
    :label: transform_hetero

    s_{ij,\text{het}}(\mu_{ij}) = s_{ij}(\nu_{ij})

The parameter :math:`a_{ij}` controls the position between atoms :math:`i` and
:math:`j` where :math:`\nu_{ij}` goes through zero, and can be used to tune
the size of the basins defined by the weight functions. Based on the covalent
bond radii, :math:`R_i`, Becke defines

.. math::
    u_{ij} = \frac{R_i-R_j}{R_i+R_j}

    a_{ij} = \frac{u_{ij}}{u_{ij}^2-1}

This choice assigns proportionally larger basins to larger atoms in the
molecule, which further improves the convergence of the numerical integrations
over the atomic grids. Note that the absolute value of :math:`a_{ij}` must be
smaller than :math:`\frac{1}{2}` to guarantee that the transform in equation
:eq:`transform_hetero` is monotonous.

In HiPart the parameter :math:`a_{ij}` is constrained to have an absolute value
smaller than 0.45 to suppress pristine behavior. The covalent radii for HiPart
are taken from [Cordero2008]_.

.. _hirshfeld:

Hirshfeld scheme
^^^^^^^^^^^^^^^^

Hirshfeld [Hirshfeld1977]_ proposed a partitioning scheme where pro-atomic
densities are derived from computations on neutral atoms by simply averaging the
atomic density over the angular degrees of freedom,

.. math:: \rho_i^{\text{pro}}(|\mathbf{r} - \mathbf{r}_i|) = \rho_i^{\text{pro}}(r) = \int d\Omega \rho_{i,N=Z}^{\text{atom}}(r,\Omega),

where :math:`\Omega` represents the angular degrees of freedom. Prior to the
application of this partitioning scheme one must setup a database of spherically
averaged atomic densities for all elements that are present in the molecule of
interest. For the sake of consistency, this needs to be carried out with the same
level of theory (and basis set) that is used for the molecular computation.

.. _hirshfeld-i:

Iterative Hirshfeld scheme
^^^^^^^^^^^^^^^^^^^^^^^^^^

The choice of neutral pro-atoms in the standard Hirshfeld scheme is somewhat
arbitrary. The Iterative Hirshfeld scheme [Bultinck2007]_ is an extension of the
original method, where one seeks for pro-atomic densities that have the same
number of electrons as the atomic partitions in the molecule.

Bultinck et al. introduce a pro-atomic function with an additional parameter,
:math:`N_i`, ie.e the fractional number of electrons in the pro-atomic density.
For integer values of this parameter, the pro-atomic density is just the
spherical average of the corresponding atom in vacuum:

.. math:: \rho_i^{\text{pro}}(r;N_i) = \int d\Omega \rho_{i,N=N_i}^{\text{atom}}(r,\Omega).

For non-integer values of the parameter N_i, the pro-atomic density is a linear
interpolation between the two `neighboring` integer-charged atoms:

.. math:: \rho_i^{\text{pro}}(r;N_i) = (\mathrm{ceil}(N_i)-N_i)\rho_i^{\text{pro}}(r;\mathrm{floor}(N_i)) +
                            (N_i-\mathrm{floor}(N_i))\rho_i^{\text{pro}}(r;\mathrm{ceil}(N_i))

The values :math:`N_i` are obtained in an iterative procedure. Initially, they
are all set to zero, and one computes the populations just like in the original
Hirshfeld scheme. In the subsequent iterations the parameters :math:`N_i` are
set to the populations from the previous iteration and one uses these pro-atoms
to compute the population for the next iteration. This is repeated until the
atomic populations converge, i.e. when the maximum absolute value of the
difference in atomic populations between two iterations drops below a predefined
threshold.

Before one can use the Iterative Hirshfeld methods, one must first construct
a database of pro-atomic densities for all the elements in the molecule under
scrutiny. For each element one must compute different charge states.


Iterative Stockholder Analysis scheme
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ISA scheme is another extension to the original Hirshfeld method where one
tries to construct spherically symatric pro-atoms that are as similar as
possible to the atomic partitions in the molecule. [Lillestolen2008]_

The initial pro-atoms are constructed by taking the minimal molecular electron
density as a function of the distance from the nucleus. For numerical reasons
this minimal value constrained to be non-zero:

.. math:: \rho_i^{\text{pro},(0)}(|\mathbf{r} - \mathbf{r}_i|) = \max(\epsilon, \min_{\Omega_i} \rho_{\text{mol}}(|\mathbf{r} - \mathbf{r}_i|,\Omega_i))

where :math:`epsilon` is a small positive number and :math:`\Omega_i` are the
angular degrees of freedom of the spherical coordinate system centered at atom
:math:`i`. In each ISA iteration :math:`k`, the new pro-atoms are taken to be
the spherical average of the atomic densities from the previous iteration.

.. math:: \rho_i^{\text{pro},(k+1)}(|\mathbf{r} - \mathbf{r}_i|) = \int d \Omega_i w_i^{(k)}(|\mathbf{r} - \mathbf{r}_i|) \rho_{\text{mol}}(|\mathbf{r} - \mathbf{r}_i|,\Omega_i)

This is again repeated until the atomic populations converge. Note that this
scheme does not depend on a database of atomic densities.


Atomic properties derived from the partitioned density
------------------------------------------------------

TODO

Charges, Dipoles & Multipoles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Gross versus net charges
^^^^^^^^^^^^^^^^^^^^^^^^

TODO

Atomic properties derived from the partitioned orbitals
-------------------------------------------------------

TODO

Spin charges
^^^^^^^^^^^^

TODO

Bond orders
^^^^^^^^^^^

TODO

Overlap elements
^^^^^^^^^^^^^^^^

TODO

Electrostatic Potential Fitting
-------------------------------

TODO
