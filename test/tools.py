# HiPart is a software toolkit to analyse molecular densities with the hirshfeld partitioning scheme.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from hipart.cache import HirshfeldICache
from hipart.context import Context
from hipart.tools import compute_stockholder_weights, load_charges, load_dipoles
from molmod import angstrom

import unittest


__all__ = ["ToolsTestCase"]


class FakeOptions(object):
    def __init__(self, density, reference, lebedev, mol_lebedev, clean, threshold, max_iter, fix_total_charge):
        self.density = density
        self.reference = reference
        self.lebedev = lebedev
        self.mol_lebedev = mol_lebedev
        self.clean = clean
        self.threshold = threshold
        self.max_iter = max_iter
        self.fix_total_charge = fix_total_charge


class ToolsTestCase(unittest.TestCase):
    def test_compute_stockholder_weights(self):
        options = FakeOptions("scf", None, 50, 50, 1, 1e-4, 500, True)
        context = Context("input/hcl.fchk", options)
        cache = HirshfeldICache.new_from_args(context, ["input/HF-STO-3G-densities.txt"])
        cache.do_atom_grids()
        cache.do_partitions()
        h0 = compute_stockholder_weights(0, cache.pro_atom_fns, context.num_lebedev, cache.atom_grid_distances, k=1)
        #print h0[:50]
        #print h0[4000:4050]
        #print h0[-50:]
        #print h0.shape
        h1 = compute_stockholder_weights(1, cache.pro_atom_fns, context.num_lebedev, cache.atom_grid_distances, k=1)
        #print h1[:50]
        #print h1[4000:4050]
        #print h1[-50:]
        #print h1.shape
        #print h1+h0
        #print cache.get_rs(0,0)[-1]/angstrom
        self.assertAlmostEqual(abs(h1+h0-1).max(), 0, 10)
        cache.do_charges()
        cache.do_dipoles()
        context.clean()

    def test_load_charges(self):
        charges = load_charges("input/hcl.hipart/hirshi_charges.txt")
        self.assertAlmostEqual(charges[0], -0.19470)
        self.assertAlmostEqual(charges[1],  0.19470)

    def test_load_dipoles(self):
        dipoles = load_dipoles("input/hcl.hipart/hirshi_dipoles.txt")
        self.assertAlmostEqual(dipoles[0,1], -0.00004)
        self.assertAlmostEqual(dipoles[1,2], -0.03344)

