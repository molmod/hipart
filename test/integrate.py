# HiPart is a tool to analyse molecular densities with the hirshfeld partitioning scheme
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


from hipart.core import *

import unittest, numpy, pylab


__all__ = ["IntegrateTestCase"]


class IntegrateTestCase(unittest.TestCase):
    def test_int1(self):
        x = numpy.arange(0.0, 1.0, 1e-4)
        y = x**2
        self.assertAlmostEqual(integrate_equi(x,y), 1.0/3,3)

    def test_gauss_potential(self):
        from scipy.special import erf
        r_low = 1e-5
        r_high = 10.0
        r_steps = 500
        r_factor = (r_high/r_low)**(1.0/(r_steps-1))
        rs = r_low*r_factor**numpy.arange(r_steps)
        beta = 0.5
        norm = (numpy.pi*beta**2)**(3.0/2.0)
        ys = numpy.exp(-(rs/beta)**2)/norm
        self.assertAlmostEqual(integrate(rs, ys*rs**2*4*numpy.pi), 1.0, 3)
        vs = erf(rs/beta)/rs
        vs_numer1 = -cumul_integrate(rs, cumul_integrate(rs, ys*4*numpy.pi*rs**2)/rs**2)

        qs = cumul_integrate(rs, ys*4*numpy.pi*rs**2)
        rs_inv = 1/rs
        vs_numer2 = cumul_integrate(rs_inv[::-1], qs[::-1])[::-1]

        pylab.clf()
        pylab.plot(rs, vs, label="analytic")
        pylab.plot(rs, vs_numer1, label="numerical 1")
        pylab.plot(rs, vs_numer2, label="numerical 2")
        pylab.legend(loc=0)
        pylab.savefig("output/gauss_potential.png")



