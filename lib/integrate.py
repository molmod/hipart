# HiPart is a software toolkit to analyse molecular densities with the hirshfeld partitioning scheme.
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


from hipart.spline import CubicSpline

import numpy


__all__ = [
    "integrate_equi", "integrate_log", "cumul_integrate_equi",
    "cumul_integrate_log", "integrate_lebedev",
]


def integrate_equi(u, g):
    """Computes the integral of g(u). (prefers equidistant grids)

    u is a one-dimensional array of grid points and g is an array of the same
    shape that contains the function values. An auxiliary cubic spline is
    generated for the integration. This means that equidistant grids are
    preferential but not mandatory.
    """
    s = CubicSpline(u, g)
    return s.integrate()

def integrate_log(r, f):
    """Computes the integral of f(r). (prefers logarithmic grids)

    r is a one-dimensional array of grid points and f is an array of the same
    shape that contains the function values. The integrand is computed in
    u = ln(r) with the function integrate_equi. This means that logarithmic
    grids are preferential but not mandatory."""
    u = numpy.log(r)
    g = f*r
    return integrate_equi(u,g)

def cumul_integrate_equi(u,g):
    """Computes the integral of g(u) from u[0] to all u[i]. (prefers equidistant grids)

    See integrate_equi for more details.
    """
    s = CubicSpline(u, g)
    return s.cumul_integrate()

def cumul_integrate_log(r, f):
    """Computes the integral of f(r) from r[0] to all r[i]. (prefers logarithmic grids)

    See integrate_log for more details.
    """
    u = numpy.log(r)
    g = f*r
    return cumul_integrate_equi(u,g)

def integrate_lebedev(weights, f):
    """Integrates f(phi,theta) where f is evaluated at lebedev grid points.

    The argument weights is a one-dimensional array that contains the lebedev
    weights. f is also a one-dimensional array that contains the evaluated
    function at the corresponding grid points. f can contain data for multiple
    spheres (just concatenated). In the latter case, an array of integrals over
    each sphere is returned.
    """
    num_lebedev = len(weights)
    f = f.reshape((-1,num_lebedev))
    return (f*weights).sum(axis=1)*4*numpy.pi


