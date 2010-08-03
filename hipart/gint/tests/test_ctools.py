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


from hipart.gint import reorder_dmat, add_orbitals_to_dmat, dmat_to_full

import numpy


def test_reorder_dmat1():
    # This test uses a very simple small matrix.
    dmat = numpy.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    permutation = numpy.array([1,2,0])
    reorder_dmat(dmat, permutation)
    expected_dmat = numpy.array([2.0, 4.0, 5.0, 1.0, 3.0, 0.0])
    assert(abs(dmat-expected_dmat).max() < 1e-10)


def test_reorder_dmat2():
    # This is a test based on a randomly generated matrix.
    num_dof = 10
    size = (num_dof*(num_dof+1))/2
    dmat = numpy.random.normal(0, 1, size)
    permutation = numpy.random.permutation(num_dof)
    square = numpy.zeros((num_dof, num_dof), float)
    j = 0
    for i in xrange(num_dof):
        square[i,0:i+1] = dmat[j:j+i+1]
        square[0:i+1,i] = dmat[j:j+i+1]
        j += i+1
    square = square[permutation].transpose()[permutation]
    expected_dmat = numpy.zeros(size, float)
    j = 0
    for i in xrange(num_dof):
        expected_dmat[j:j+i+1] = square[i,0:i+1]
        j += i+1
    reorder_dmat(dmat, permutation)
    assert(abs(dmat-expected_dmat).max() < 1e-10)


def test_add_orbitals_to_dmat1():
    orbitals = numpy.array([[1, 2], [3, 4]])
    expected_dmat = numpy.array([10, 14, 20])
    dmat = numpy.zeros(3, float)
    add_orbitals_to_dmat(orbitals, dmat)
    assert(abs(dmat - expected_dmat).max() < 1e-10)


def test_add_orbitals_to_dmat2():
    num_dof = 10
    num_orbitals = 5
    size = (num_dof*(num_dof+1))/2
    orbitals = numpy.random.normal(0, 1, (num_orbitals, num_dof))
    dmat = numpy.zeros(size, float)
    add_orbitals_to_dmat(orbitals, dmat)
    expected_dmat = numpy.dot(orbitals.transpose(), orbitals)
    expected_dmat = numpy.concatenate([expected_dmat[i,:i+1] for i in xrange(num_dof)])
    assert(abs(dmat - expected_dmat).max() < 1e-10)


def test_dmat_to_full1():
    dmat = numpy.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    full = numpy.zeros((3,3), float)
    dmat_to_full(dmat, full)
    expected_full = numpy.array([[1.0, 2.0, 4.0], [2.0, 3.0, 5.0], [4.0, 5.0, 6.0]])
    error = abs(full - expected_full).max()
    assert(error < 1e-10)

def test_dmat_to_full2():
    num_dof = 10
    size = (num_dof*(num_dof+1))/2
    dmat = numpy.random.normal(0, 1, size)
    full = numpy.zeros((num_dof,num_dof), float)
    dmat_to_full(dmat, full)
    expected_full = numpy.zeros((num_dof,num_dof), float)
    counter = 0
    for i in xrange(num_dof):
        row = dmat[counter:counter+i+1]
        expected_full[i,:i+1] = row
        expected_full[:i+1,i] = row
        counter += i+1
    error = abs(full - expected_full).max()
    assert(error < 1e-10)
