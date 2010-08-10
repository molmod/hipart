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


from hipart.csext import spline_construct, spline_int
from hipart.lebedev_laikov import get_grid as get_lebedev_grid

from molmod import Rotation

import numpy, os, glob


__all__ = [
    "Grid", "AtomicGrid",
    "RBaseIntGrid", "RLogIntGrid", "REquiIntGrid", "get_rgrid_from_description",
    "ABaseIntGrid", "ALebedevIntGrid"
]


class Grid(object):
    def __init__(self, prefix, work, points, dump=True):
        # check arguments
        if len(points.shape) != 2 or points.shape[1] != 3:
            raise TypeError("The points argument must be an array with 3D coordinates, i.e. shape (n,3)")
        self.prefix = prefix
        self.points = points
        self.work = work
        if dump:
            work.dump(prefix, points)

    size = property(lambda self: len(self.points))

    @classmethod
    def from_prefix(cls, prefix, work):
        points = work.load(prefix, (-1,3))
        if points is not None:
            return cls(prefix, work, points, False)
        else:
            return None

    def load(self, suffix, shape=None):
        return self.work.load("%s_%s" % (self.prefix, suffix), shape)

    def dump(self, suffix, array, ignore=False):
        self.work.dump("%s_%s" % (self.prefix, suffix), array, ignore)


class AtomicGrid(Grid):
    @classmethod
    def from_parameters(cls, prefix, work, center, rgrid, agrid):
        points = agrid.generate_points(center, rgrid.rs)
        return cls(prefix, work, points)


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


class ABaseIntGrid(object):
    def generate_points(self, center, rs):
        raise NotImplementedError

    def integrate(self, integrand):
        raise NotImplementedError

    def minimum(self, function):
        raise NotImplementedError

    def get_description(self):
        raise NotImplementedError


class ALebedevIntGrid(ABaseIntGrid):
    def __init__(self, num_lebedev, do_random=True):
        self.num_lebedev = num_lebedev
        self.do_random = do_random
        self.lebedev_xyz, self.lebedev_weights = get_lebedev_grid(num_lebedev)

    def generate_points(self, center, rs):
        points = numpy.zeros((self.num_lebedev*len(rs),3), float)
        start = 0
        for r in rs:
            end = start + self.num_lebedev
            xyz = self.lebedev_xyz
            if self.do_random:
                rot = Rotation.random()
                xyz = numpy.dot(xyz,rot.r)
            points[start:end] = r*xyz+center
            start = end
        return points

    def integrate(self, integrand):
        integrand = integrand.reshape((-1,self.num_lebedev))
        return (integrand*self.lebedev_weights).sum(axis=1)*(4*numpy.pi)

    def minimum(self, function):
        return function.reshape((-1,self.num_lebedev)).min(axis=1)

    def get_description(self):
        return "ALebedevIntGrid(%i)" % self.num_lebedev
