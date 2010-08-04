
HiPart documentation
====================

Introduction
------------

HiPart is a program to analyze the electronic structure of molecules with
fuzzy-atom partitioning methods. It now supports four schemes to define atomic
partitions: the Becke scheme [Becke1988]_, The Hirshfeld scheme
[Hirshfeld1977]_, the Iterative Hirshfeld scheme [Bultinck2007]_, and the
Iterative Stockholder Analysis [Lillestolen2008]_.

Within each scheme the following quantities can be computed: atomic charge,
atomic dipoles, quality of charges and dipoles with respect to the ESP, the
atomic multipole expansion, net and overlap populations, bond orders, spin
charges, atomic overlap matrices in the orbital basis and in the basis of
contracted Gaussians.

HiPart is developed at the Center for Molecular Modeling since October 2008 by
Toon Verstraelen. It was originally conceived as an autodidactic project to get
familiar with fuzzy atom partitioning in general, and it was used to generate
Iterative Hirshfeld charges for a benchmark study of the EEM and SQE model.
[Verstraelen2009]_. Over the past two years, the program evolved into a more
efficient and more user-friendly toolkit. As of August 2010 HiPart is publicly
available as a Open Source software. To date their is no official stable release
of HiPart and it is recommend to upgrade regularly to the latest development
release.

HiPart is open source software distributed under the conditions of the `GPL v3
license <http://molmod.ugent.be/code/wiki/GPL_License_v3>`_.


Tutorial
--------

.. toctree::
   :maxdepth: 3

   tut_install
   tut_theory
   tut_usage


Library reference
-----------------

The library reference is automatically generated from the source code. As the
docstrings in the source code get more and more complete, this part will become
more and more useful.

.. toctree::
   :maxdepth: 3

   ref_hipart
   ref_gint


References
----------

.. toctree::
   :maxdepth: 3

   references


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
