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

import numpy


def test_blind1():
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

def test_blind2():
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
