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
from hipart.log import log

import numpy


__all__ = [
    "LinearSpline", "CubicSpline",
    "RBaseIntGrid", "RLogIntGrid", "REquiIntGrid", "get_rgrid_from_description",
]


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


class RBaseIntGrid(object):
    def __init__(self, rs):
        self.rs = rs
        self._weights_map = {}

    def _get_weights_aux(self, us, vs):
        ys = numpy.zeros(len(us), float)
        ds = numpy.zeros(len(us), float)
        ws = numpy.zeros(len(us), float)
        for i in xrange(len(us)):
            ys[:] = 0
            if vs is None:
                ys[i] = 1
            else:
                ys[i] = vs[i]
            spline_construct(us,ys,ds)
            ws[i] = spline_int(us,ys,ds)
        return ws

    def _get_weights_low(self, size):
        raise NotImplementedError

    def get_weights(self, size=None):
        if size is None:
            size = len(self.rs)
        weights = self._weights_map.get(size)
        if weights is None:
            weights = self._get_weights_low(size)
            self._weights_map[size] = weights
        return weights

    def get_description(self):
        raise NotImplementedError

    def integrate(self, integrand):
        w = self.get_weights(len(integrand))
        return numpy.dot(w, integrand)

    def integrate_cumul(self, integrand):
        # WARNING: this is very inaccurate, do not use for delicate computations
        w = self.get_weights(len(integrand))
        return numpy.cumsum(w*integrand)



class REquiIntGrid(RBaseIntGrid):
    def __init__(self, r_low, r_high, steps):
        if r_low >= r_high:
            raise ValueError("The argument r_high must be strictly larger than r_low")
        if steps <= 0:
            raise ValueError("The argument steps must be strictly positive")
        self.r_low = r_low
        self.r_high = r_high
        self.steps = steps
        delta = (r_high - r_low)/(steps - 1)
        rs = delta*numpy.arange(0,steps) + r_low
        RBaseIntGrid.__init__(self, rs)

    def _get_weights_low(self, size):
        return self._get_weights_aux(self.rs[:size], None)

    def get_description(self):
        return "REquiIntGrid(%#.16e,%#.16e,%i)" % (self.r_low,self.r_high,self.steps)


class RLogIntGrid(RBaseIntGrid):
    def __init__(self, r_low, r_high, steps):
        if r_low <= 0:
            raise ValueError("The argument r_low must be strictly positive")
        if r_high <= 0:
            raise ValueError("The argument r_high must be strictly positive")
        if r_low >= r_high:
            raise ValueError("The argument r_high must be strictly larger than r_low")
        if steps <= 0:
            raise ValueError("The argument steps must be strictly positive")
        self.r_low = r_low
        self.r_high = r_high
        self.steps = steps
        ratio = (r_high/r_low)**(1.0/(steps-1))
        alpha = numpy.log(ratio)
        rs = r_low*numpy.exp(alpha*numpy.arange(0,steps))
        RBaseIntGrid.__init__(self, rs)

    def _get_weights_low(self, size):
        rs = self.rs[:size]
        us = numpy.log(rs)
        vs = rs
        return self._get_weights_aux(us, vs)

    def get_description(self):
        return "RLogIntGrid(%#.16e,%#.16e,%i)" % (self.r_low,self.r_high,self.steps)


def get_rgrid_from_description(s):
    pos = s.find("(")
    name = s[:pos]
    arg_strs = s[pos+1:-1].split(",")
    args = []
    for arg in arg_strs:
        if "." in arg:
            args.append(float(arg))
        else:
            args.append(int(arg))
    Classes = {
        "REquiIntGrid": REquiIntGrid,
        "RLogIntGrid": RLogIntGrid,
    }
    return Classes[name](*args)
