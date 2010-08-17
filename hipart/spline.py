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


from hipart.ext import spline_construct, spline_eval
from hipart.log import log

import numpy


__all__ = ["LinearSpline", "CubicSpline"]


class LinearSpline(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
        a = (y[1:] - y[:-1])/(x[1:] - x[:-1])
        self.a = numpy.concatenate(([0], a, [0]))
        self.r0 = numpy.concatenate(([0], x))
        self.b = numpy.concatenate(([y[0]], y))

    def __call__(self, x):
        # simple linear interpolation
        i = self.x.searchsorted(x)
        return (x-self.r0[i])*self.a[i]+self.b[i]


class CubicSpline(object):
    def __init__(self, x, y):
        if len(x) != len(y):
            raise ValueError("x and y must have the same length.")
        if ((x[1:] - x[:-1]) < 0).any():
            raise ValueError("x must be a sorted array.")
        self.x = numpy.array(x, float, copy=False)
        self.y = numpy.array(y, float, copy=False)
        self.d = numpy.zeros(len(self.x), float)
        spline_construct(self.x,self.y,self.d)

    def __call__(self, x):
        y = numpy.zeros(len(x), float)
        spline_eval(self.x,self.y,self.d,x,y)
        return y
