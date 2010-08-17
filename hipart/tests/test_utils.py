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


from hipart.ext import grid_distances

import numpy


def test_grid_distances():
    N = 10
    p = numpy.random.normal(0,1,(N,3))
    c = numpy.random.normal(0,1,(3,))
    d1 = numpy.sqrt(((p - c)**2).sum(axis=1))
    d2 = numpy.zeros(d1.shape, float)
    grid_distances(c, p, d2)
    error = abs(d1 - d2).max()
    assert(error < 1e-5)
