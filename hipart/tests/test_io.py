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


from hipart.io import *

import tempfile, shutil, os, numpy


hirshi_charges_txt="""\
number of atoms: 2
  i        Z      Charge
-----------------------------
  1   F    9     -0.19457
  2   H    1      0.19457
"""

def test_load_atom_scalars():
    tmpdir = tempfile.mkdtemp("hipart-load-charges")
    fn_txt = os.path.join(tmpdir, "hirshi_charges.txt")
    f = open(fn_txt, "w")
    f.write(hirshi_charges_txt)
    f.close()
    charges = load_atom_scalars(fn_txt)
    assert(charges.shape == (2,))
    assert(abs(charges[0] - (-0.19438)) < 1e-3)
    assert(abs(charges[1] - 0.19438) < 1e-3)
    shutil.rmtree(tmpdir)


def test_dump_atom_scalars():
    tmpdir = tempfile.mkdtemp("hipart-dump-charges")
    fn_txt = os.path.join(tmpdir, "foo_charges.txt")
    charges = [-0.3, 0.5]
    numbers = [2, 3]
    dump_atom_scalars(fn_txt, charges, "Scalar", numbers)
    check = load_atom_scalars(fn_txt)
    assert(len(charges)==len(check))
    assert(abs(charges-check).max() < 1e-5)
    dump_atom_scalars(fn_txt, charges, "Scalar")
    check = load_atom_scalars(fn_txt)
    assert(len(charges)==len(check))
    assert(abs(charges-check).max() < 1e-5)
    shutil.rmtree(tmpdir)


hirshi_dipoles_txt="""\
number of atoms: 2
  i        Z     Dipole-X     Dipole-Y     Dipole-Z      Dipole
------------------------------------------------------------------
  1   F    9      0.00000      0.00000     -0.06956      0.06956
  2   H    1      0.00003     -0.00001     -0.03371      0.03371
"""

def test_load_atom_vectors():
    tmpdir = tempfile.mkdtemp("hipart-load-dipoles")
    fn_txt = os.path.join(tmpdir, "hirshi_dipoles.txt")
    f = open(fn_txt, "w")
    f.write(hirshi_dipoles_txt)
    f.close()
    dipoles = load_atom_vectors(fn_txt)
    assert(dipoles.shape == (2,3))
    assert(abs(dipoles[0,1] -  0.00002) < 1e-3)
    assert(abs(dipoles[1,2] - (-0.03385)) < 1e-3)
    shutil.rmtree(tmpdir)


def test_dump_atom_vectors():
    tmpdir = tempfile.mkdtemp("hipart-dump-dipoles")
    fn_txt = os.path.join(tmpdir, "hirshi_dipoles.txt")
    dipoles = numpy.array([
        [1.0, 2.0, 3.0],
        [-1.0, 1.0, 0.0],
    ])
    numbers = [1, 6]
    dump_atom_vectors(fn_txt, dipoles, "Vector", numbers)
    check = load_atom_vectors(fn_txt)
    assert(dipoles.shape==check.shape)
    assert(abs(dipoles-check).max() < 1e-5)
    dump_atom_vectors(fn_txt, dipoles, "Vector")
    check = load_atom_vectors(fn_txt)
    assert(dipoles.shape==check.shape)
    assert(abs(dipoles-check).max() < 1e-5)
    shutil.rmtree(tmpdir)


hirshi_bond_orders_txt = """\
number of atoms: 8
   Bond order   |     1  C       2  O       3  O       4  H       5  C       6  H       7  H       8  H
----------------+-----------------------------------------------------------------------------------------
  1   C    6    |    0.00000    1.19387    1.82032    0.04743    0.97009    0.05155    0.05472    0.05476
  2   O    8    |    1.19387    0.00000    0.33499    0.78470    0.22559    0.01126    0.02059    0.02064
  3   O    8    |    1.82032    0.33499    0.00000    0.04033    0.25623    0.02708    0.01639    0.01638
  4   H    1    |    0.04743    0.78470    0.04033    0.00000    0.01245    0.00118    0.00084    0.00085
  5   C    6    |    0.97009    0.22559    0.25623    0.01245    0.00000    0.94421    0.94237    0.94216
  6   H    1    |    0.05155    0.01126    0.02708    0.00118    0.94421    0.00000    0.05648    0.05619
  7   H    1    |    0.05472    0.02059    0.01639    0.00084    0.94237    0.05648    0.00000    0.05849
  8   H    1    |    0.05476    0.02064    0.01638    0.00085    0.94216    0.05619    0.05849    0.00000
"""

def test_load_atom_matrix():
    tmpdir = tempfile.mkdtemp("hipart-load-matrix")
    fn_txt = os.path.join(tmpdir, "hirshi_bond_orders.txt")
    f = open(fn_txt, "w")
    f.write(hirshi_bond_orders_txt)
    f.close()
    matrix = load_atom_matrix(fn_txt)
    assert(matrix.shape == (8,8))
    assert(abs(matrix[0,0]) < 1e-4)
    assert(abs(matrix[1,2] - 0.33499) < 1e-4)
    assert(abs(matrix[-1,2] - 0.01638) < 1e-4)
    shutil.rmtree(tmpdir)


def test_dump_atom_matrix():
    tmpdir = tempfile.mkdtemp("hipart-dump-matrix")
    fn_txt = os.path.join(tmpdir, "hirshi_matrix.txt")
    matrix = numpy.random.normal(0,1,(5,5))
    numbers = numpy.random.randint(1,10,5)
    dump_atom_matrix(fn_txt, matrix, "Matrix", numbers)
    check = load_atom_matrix(fn_txt)
    assert(matrix.shape==check.shape)
    assert(abs(matrix-check).max() < 1e-5)
    dump_atom_matrix(fn_txt, matrix, "Matrix")
    check = load_atom_matrix(fn_txt)
    assert(matrix.shape==check.shape)
    assert(abs(matrix-check).max() < 1e-5)
    shutil.rmtree(tmpdir)


some_fields_txt = """\
number of atoms: 2
number of labels: 3
   bla-bla-bla  |     one        two        three
----------------+----------------------------------
  1   C    6    |    1.50000    3.00000    3.00000
  2   O    8    |    4.00000    2.00000    5.00000
"""

def test_load_atom_fields():
    tmpdir = tempfile.mkdtemp("hipart-load-fields")
    fn_txt = os.path.join(tmpdir, "some_fields.txt")
    f = open(fn_txt, "w")
    f.write(some_fields_txt)
    f.close()
    table, labels = load_atom_fields(fn_txt)
    assert(labels == ["one", "two", "three"])
    assert(table.shape == (2,3))
    assert(abs(table[0,0] - 1.5) < 1e-4)
    assert(abs(table[1,1] - 2.0) < 1e-4)
    assert(abs(table[-1,2] - 5.0) < 1e-4)
    shutil.rmtree(tmpdir)


def test_dump_atom_fields():
    tmpdir = tempfile.mkdtemp("hipart-dump-fields")
    fn_txt = os.path.join(tmpdir, "fields.txt")
    table = numpy.random.normal(0,1,(5,8))
    numbers = numpy.random.randint(1,10,5)
    labels = list("abcdefgh")
    dump_atom_fields(fn_txt, table, labels, "Table", numbers)
    check_table, check_labels = load_atom_fields(fn_txt)
    assert(labels==check_labels)
    assert(table.shape==check_table.shape)
    assert(abs(table-check_table).max() < 1e-5)
    dump_atom_fields(fn_txt, table, labels, "Table")
    check_table, check_labels = load_atom_fields(fn_txt)
    assert(labels==check_labels)
    assert(table.shape==check_table.shape)
    assert(abs(table-check_table).max() < 1e-5)
    shutil.rmtree(tmpdir)
