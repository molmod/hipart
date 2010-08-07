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


from hipart.grids import *

import numpy


def test_radial_integration_equi():
    grid = REquiIntGrid(0.0,1.0,100)
    x = grid.rs
    weights = grid.get_weights(len(x))
    #
    y = numpy.sin(x)
    int1 = numpy.dot(y,weights)
    int2 = 1 - numpy.cos(x[-1])
    assert(abs(int1-int2)<1e-7)
    #
    y = x*x
    int1 = numpy.dot(y,weights)
    int2 = x[-1]**3/3.0
    assert(abs(int1-int2)<1e-7)
    #
    y = numpy.exp(-x)
    int1 = numpy.dot(y,weights)
    int2 = 1.0 - numpy.exp(-x[-1])
    assert(abs(int1-int2)<1e-7)

def test_radial_integration_log():
    grid = RLogIntGrid(1e-5,1e2,100)
    x = grid.rs
    weights = grid.get_weights(len(x))
    #
    y = x*numpy.exp(-x)
    int1 = numpy.dot(y,weights)
    int2 = 1.0
    error = abs(int1-int2)
    assert(error<1e-10)
    #
    y = x*x*x*numpy.exp(-x)
    int1 = numpy.dot(y,weights)
    int2 = 6.0
    error = abs(int1-int2)
    assert(error<1e-14)
