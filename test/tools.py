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






from hipart.cache import HirshfeldCache
from hipart.context import Context
from hipart.tools import load_charges, load_dipoles, dump_charges
from molmod import angstrom

import unittest


__all__ = ["ToolsTestCase"]


class FakeOptions(object):
    def __init__(self, lebedev, mol_lebedev, clean, threshold, max_iter, fix_total_charge):
        self.lebedev = lebedev
        self.mol_lebedev = mol_lebedev
        self.clean = clean
        self.threshold = threshold
        self.max_iter = max_iter
        self.fix_total_charge = fix_total_charge


class ToolsTestCase(unittest.TestCase):
    def test_compute_stockholder_weights(self):
        options = FakeOptions(110, 50, 1, 1e-4, 500, True)
        context = Context("input/hcl.fchk", options)
        cache = HirshfeldCache.new_from_args(context, ["input/HF-STO-3G-densities.txt"])
        cache.do_atgrids()
        cache.do_proatomfns()
        h0 = cache._compute_atweights(cache.atgrids[0], 0)
        #print h0[:50]
        #print h0[4000:4050]
        #print h0[-50:]
        #print h0.shape
        h1 = cache._compute_atweights(cache.atgrids[0], 1)
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
        self.assertAlmostEqual(charges[0], -0.19438, 3)
        self.assertAlmostEqual(charges[1],  0.19438, 3)

    def test_dump_charges(self):
        fn_txt = "output/foo_charges.txt"
        charges = [-0.3, 0.5]
        numbers = [2, 3]
        dump_charges(fn_txt, charges, numbers)
        check = load_charges(fn_txt)
        self.assertEqual(len(charges), len(check))
        self.assertAlmostEqual(charges[0], check[0])
        self.assertAlmostEqual(charges[1], check[1])
        dump_charges(fn_txt, charges)
        check = load_charges(fn_txt)
        self.assertEqual(len(charges), len(check))
        self.assertAlmostEqual(charges[0], check[0])
        self.assertAlmostEqual(charges[1], check[1])

    def test_load_dipoles(self):
        dipoles = load_dipoles("input/hcl.hipart/hirshi_dipoles.txt")
        self.assertAlmostEqual(dipoles[0,1],  0.00002, 3)
        self.assertAlmostEqual(dipoles[1,2], -0.03385, 3)


