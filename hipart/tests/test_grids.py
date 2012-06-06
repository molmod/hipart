# -*- coding: utf-8 -*-
# HiPart is a program to analyze the electronic structure of molecules with
# fuzzy-atom partitioning methods.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
#--


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


def test_rgrid_from_description():
    grid1 = RLogIntGrid(3.7794522678425048e-05,3.7794522678425039e+01,100)
    grid2 = get_rgrid_from_description("RLogIntGrid(3.7794522678425048e-05,3.7794522678425039e+01,100)")
    error = abs(grid1.rs - grid2.rs).max()
    assert(error<1e-12)

def test_rgrid_sanity():
    grid = REquiIntGrid(0.0,1.0,100)
    assert(grid.rs[0] == 0.0)
    assert(grid.rs[-1] == 1.0)
    assert(len(grid.rs) == 100)
    grid = RLogIntGrid(1e-5,1e2,100)
    assert(grid.rs[0] == 1e-5)
    assert(abs(grid.rs[-1] - 1e2) < 1e-10)
    assert(len(grid.rs) == 100)
