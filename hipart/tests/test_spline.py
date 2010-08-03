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


from hipart.csext import *
from hipart.spline import *

import numpy


def test_blind():
    # see whether we get segfaults or not.
    x = numpy.arange(0,1,0.1)
    y = numpy.random.normal(0,1,len(x))
    d = numpy.zeros(len(x), float)
    spline_construct(x,y,d)
    x_new = numpy.random.uniform(-1, 2, 50)
    y_new = numpy.zeros(len(x_new), float)
    spline_eval(x,y,d,x_new,y_new)
    yint = numpy.zeros(len(y), float)
    spline_cumul_int(x,y,d,yint)

def test_plot1():
    x = numpy.arange(0, 1.005, 0.1)*numpy.pi*2
    y = numpy.sin(x)
    yint = -numpy.cos(x)
    yint -= yint[0]
    d = numpy.zeros(len(x), float)
    spline_construct(x,y,d)
    x_new = numpy.arange(0,1.0005, 0.01)*numpy.pi*2
    y_new = numpy.zeros(len(x_new), float)
    spline_eval(x,y,d,x_new,y_new)
    yint_spline = numpy.zeros(len(y), float)
    spline_cumul_int(x,y,d,yint_spline)

def test_radial_integration_equi():
    x = numpy.arange(0.0,1.0001,0.01)
    weights = get_radial_weights_equi(x)
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
    low = 1e-5
    high = 1e2
    steps = 100
    alpha = (high/low)**(1.0/(steps-1))
    x = low*alpha**numpy.arange(steps)
    weights = get_radial_weights_log(x)
    #print x
    #print weights
    #
    y = x*numpy.exp(-x)
    int1 = numpy.dot(y,weights)
    int2 = 1.0
    assert(abs(int1-int2)<1e-10)
    #
    y = x*x*x*numpy.exp(-x)
    int1 = numpy.dot(y,weights)
    int2 = 6.0
    assert(abs(int1-int2)<1e-14)
