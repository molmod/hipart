# HiPart is a software toolkit to analyse molecular densities with the hirshfeld partitioning scheme.
# Copyright (C) 2007 - 2010 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of HiPart.
#
# HiPart is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HiPart is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


from hipart.wavefn import *
from hipart.gint import dmat_to_full

import numpy


def test_num_filled():
    assert(get_num_filled(numpy.array([1.0, 1.0, 1.0]))==3)
    assert(get_num_filled(numpy.array([1.0, 1.0, 0.0]))==2)
    assert(get_num_filled(numpy.array([1.0, 1.0, 1e-6]))==2)


def test_compute_naturals():
    num_dof = 10
    # generate random symmetric matric
    A = numpy.random.normal(0,1,(num_dof, num_dof))
    A = A + A.transpose()
    # diagonalize to get the fake orbitals
    check_orbitals = numpy.linalg.eigh(A)[1].transpose()
    # invent occupation numbers
    check_occupations = numpy.random.normal(0, 1, num_dof)
    check_occupations.sort()
    check_occupations = check_occupations[::-1]
    # reconstruct the full density matrix
    check_foo = check_orbitals.transpose()
    full = numpy.dot(check_occupations*check_foo, check_foo.transpose())

    # intermediate test
    occupations, foo = numpy.linalg.eigh(full)
    occupations = occupations[::-1]
    foo = foo[:,::-1]
    orbitals = foo.transpose()
    error = abs(occupations - check_occupations).max()
    assert(error<1e-10)
    overlap = abs(numpy.dot(orbitals, check_orbitals.transpose()).round())
    error = abs(overlap - numpy.identity(num_dof)).max()
    assert(error<1e-10)

    # make it compact
    size = (num_dof*(num_dof+1))/2
    dmat = numpy.zeros(size, float)
    counter = 0
    for i in xrange(num_dof):
        dmat[counter:counter+i+1] = full[i,:i+1]
        counter += i+1
    check_full = numpy.zeros((num_dof, num_dof), float)
    dmat_to_full(dmat, check_full)
    error = abs(full - check_full).max()
    assert(error<1e-10)
    # run the routine to be tested
    occupations, orbitals = compute_naturals(dmat, num_dof)
    error = abs(occupations - check_occupations).max()
    assert(error<1e-10)
    overlap = abs(numpy.dot(orbitals, check_orbitals.transpose()).round())
    error = abs(overlap - numpy.identity(num_dof)).max()
    assert(error<1e-10)
