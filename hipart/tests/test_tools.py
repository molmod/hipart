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


from hipart.tools import load_charges, load_dipoles, dump_charges

from molmod import angstrom
import tempfile, shutil, os


hirshi_charges_txt="""\
number of atoms: 2
  i        Z      Charge
-----------------------------
  1   F    9     -0.19457
  2   H    1      0.19457
-----------------------------
   Q SUM          0.00000
   Q RMS          0.19457
 ESP RMS          7.88796e-03
ESP RMSD          1.75095e-03"""

def test_load_charges():
    tmpdir = tempfile.mkdtemp("hipart-load-charges")
    fn_txt = os.path.join(tmpdir, "hirshi_charges.txt")
    f = open(fn_txt, "w")
    f.write(hirshi_charges_txt)
    f.close()
    charges = load_charges(fn_txt)
    assert(abs(charges[0] - (-0.19438)) < 1e-3)
    assert(abs(charges[1] - 0.19438) < 1e-3)
    shutil.rmtree(tmpdir)


def test_dump_charges():
    tmpdir = tempfile.mkdtemp("hipart-dump-charges")
    fn_txt = os.path.join(tmpdir, "foo_charges.txt")
    charges = [-0.3, 0.5]
    numbers = [2, 3]
    dump_charges(fn_txt, charges, numbers)
    check = load_charges(fn_txt)
    assert(len(charges)==len(check))
    assert(abs(charges-check).max() < 1e-5)
    dump_charges(fn_txt, charges)
    check = load_charges(fn_txt)
    assert(len(charges)==len(check))
    assert(abs(charges-check).max() < 1e-5)
    shutil.rmtree(tmpdir)


hirshi_dipoles_txt="""\
number of atoms: 2
  i        Z     Dipole-X     Dipole-Y     Dipole-Z      Dipole
------------------------------------------------------------------
  1   F    9      0.00000      0.00000     -0.06956      0.06956
  2   H    1      0.00003     -0.00001     -0.03371      0.03371
------------------------------------------------------------------
Molecular dipole due to ...
charges (q)       0.00000      0.00000     -0.37062      0.37062
dipoles (p)       0.00003     -0.00000     -0.10327      0.10327
q and p           0.00003     -0.00000     -0.47389      0.47389
total density     0.00000      0.00000     -0.47390      0.47390
------------------------------------------------------------------
Reproduction of the external molecular ESP ...
                     RMSD             RMS       CORRELATION
charges (q)       1.75095e-03     6.22261e-03       1.00
dipoles (p)       6.14792e-03     1.78334e-03       0.98
q and p           8.26040e-04     7.99629e-03       0.99
total density                     7.88796e-03"""

def test_load_dipoles():
    tmpdir = tempfile.mkdtemp("hipart-load-dipoles")
    fn_txt = os.path.join(tmpdir, "hirshi_dipoles.txt")
    f = open(fn_txt, "w")
    f.write(hirshi_dipoles_txt)
    f.close()
    dipoles = load_dipoles(fn_txt)
    assert(abs(dipoles[0,1] -  0.00002) < 1e-3)
    assert(abs(dipoles[1,2] - (-0.03385)) < 1e-3)
    shutil.rmtree(tmpdir)
