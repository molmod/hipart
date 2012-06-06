:mod:`hipart.gint` --- (Towards a) Gaussian integrals library
=============================================================

Gint is a part of HiPart used to evaluate functions and integrals of Gaussian
functions. Gint is for the moment very limited in functionality and just does
what is needed to make HiPart work.

Gint is also an experiment to see if it is possible to let a computer program
write efficient implementations of the Gaussian integrals. (These routines are
not so funny to write by hand after all.) All C code is written by `writer`
scripts that use `SymPy <http://www.sympy.org>`_ as the computer algebra layer.
Although there are still opportunities to further optmize the generated code,
the current result is already very useful.

:mod:`hipart.gint.basis` --- Basis set representation
-----------------------------------------------------

.. automodule:: hipart.gint.basis
   :members:

One-center routines
-------------------

TODO

Two-center routines
-------------------

TODO

Auxiliary function
------------------

TODO

Tools
-----

TODO

Writer
------

TODO
