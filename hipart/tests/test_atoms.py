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


from utils import setup_hf_sto3g_gaussian

from hipart.atoms import AtomTable

import shutil, numpy


def test_atom_table():
    tmpdir, fn_fchk, fn_densities = setup_hf_sto3g_gaussian()
    atom_table = AtomTable(fn_densities)
    check_population(atom_table, 1, 0.0)
    check_population(atom_table, 1, 0.5)
    check_population(atom_table, 1, 1.0)
    check_population(atom_table, 1, 2.5)
    check_population(atom_table, 1, 5.0)
    check_population(atom_table, 8, 5.0)
    check_population(atom_table, 8, 7.5)
    check_population(atom_table, 8, 8.0)
    check_population(atom_table, 8, 9.5)
    check_population(atom_table, 8, 13.0)
    check_population(atom_table, 9, 5.0)
    check_population(atom_table, 9, 7.5)
    check_population(atom_table, 9, 8.0)
    check_population(atom_table, 9, 9.5)
    check_population(atom_table, 9, 13.0)
    check_identity(atom_table, 1, None)
    check_identity(atom_table, 1, 1)
    check_identity(atom_table, 1, 2)
    check_identity(atom_table, 8, None)
    check_identity(atom_table, 8, 6)
    check_identity(atom_table, 8, 7)
    check_identity(atom_table, 8, 8)
    check_identity(atom_table, 8, 9)
    check_identity(atom_table, 8, 10)
    check_identity(atom_table, 9, None)
    check_identity(atom_table, 9, 7)
    check_identity(atom_table, 9, 8)
    check_identity(atom_table, 9, 9)
    check_identity(atom_table, 9, 10)
    shutil.rmtree(tmpdir)

def check_population(atom_table, number, population):
    atom_profile = atom_table.records[number]
    spline = atom_profile.get_atom_fn(population)
    rs = spline.x
    rhos = spline.y
    rgrid = atom_table.rgrid
    check_population = rgrid.integrate(4*numpy.pi*rs*rs*rhos)
    error = abs(check_population - population)
    assert(error < 1e-5)

def check_identity(atom_table, number, population):
    atom_profile = atom_table.records[number]
    spline = atom_profile.get_atom_fn(population)
    if population is None:
        population = atom_profile.number
    assert(spline.y is atom_profile.records[population])
