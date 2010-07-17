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


from hipart.lebedev_laikov import grid_fns, get_grid

import numpy


def test_get_grid():
    for number in grid_fns:
        xyz, w = get_grid(number)
        assert(len(xyz.shape)==2)
        assert(len(w.shape)==1)
        assert(len(xyz)==len(w))
        assert(xyz.shape[1]==3)
        assert(abs(numpy.dot(xyz[:,0], w)) < 1e-10)
        assert(abs(numpy.dot(xyz[:,1], w)) < 1e-10)
        assert(abs(numpy.dot(xyz[:,2], w)) < 1e-10)
