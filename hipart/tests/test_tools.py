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


from hipart.tools import load_atom_scalars, load_atom_vectors, \
    dump_atom_scalars, dump_atom_vectors

from molmod import angstrom
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
    assert(abs(charges[0] - (-0.19438)) < 1e-3)
    assert(abs(charges[1] - 0.19438) < 1e-3)
    shutil.rmtree(tmpdir)


def test_dump_atom_scalars():
    tmpdir = tempfile.mkdtemp("hipart-dump-charges")
    fn_txt = os.path.join(tmpdir, "foo_charges.txt")
    charges = [-0.3, 0.5]
    numbers = [2, 3]
    dump_atom_scalars(fn_txt, charges, numbers)
    check = load_atom_scalars(fn_txt)
    assert(len(charges)==len(check))
    assert(abs(charges-check).max() < 1e-5)
    dump_atom_scalars(fn_txt, charges)
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
    dump_atom_vectors(fn_txt, dipoles, numbers)
    check = load_atom_vectors(fn_txt)
    assert(len(dipoles)==len(check))
    assert(abs(dipoles-check).max() < 1e-5)
    dump_atom_vectors(fn_txt, dipoles)
    check = load_atom_vectors(fn_txt)
    assert(len(dipoles)==len(check))
    assert(abs(dipoles-check).max() < 1e-5)
    shutil.rmtree(tmpdir)
