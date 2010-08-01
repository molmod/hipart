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


from sympy import simplify, Symbol, cos, sin, exp, Ylm, I, sqrt, pi, latex


__all__ = ["get_solid_harmonics"]


def get_solid_harmonics(shell, xyz):
    x, y, z = xyz
    theta = Symbol("theta", real=True)
    phi = Symbol("phi", real=True)
    r = Symbol("r", real=True)
    result = []
    for m in range(shell+1):
        if m > 0:
            part_plus = (Ylm(shell,m,theta,phi)*r**shell).expand()
            part_plus = part_plus.subs(exp(I*phi),(x+I*y)/sin(theta)/r)
            part_min = (Ylm(shell,-m,theta,phi)*r**shell).expand()
            part_min = part_min.subs(exp(-I*phi),(x-I*y)/sin(theta)/r)
            sym = simplify(((-1)**m*part_plus+part_min)/sqrt(2))
            sym = sym.subs(r*cos(theta),z)
            sym = sym.subs(r**2,x**2+y**2+z**2)
            sym = sym.subs(cos(theta)**2,1-sin(theta)**2)
            sym = simplify(sym)
            asym = simplify(((-1)**m*part_plus-part_min)/I/sqrt(2))
            asym = asym.subs(r*cos(theta),z)
            asym = asym.subs(r**2,x**2+y**2+z**2)
            asym = asym.subs(cos(theta)**2,1-sin(theta)**2)
            asym = simplify(asym)
            result.append(sym)
            result.append(asym)
        else:
            part = (Ylm(shell,0,theta,phi)*r**shell).expand()
            part = part.subs(r*cos(theta),z)
            part = part.subs(r**2,x**2+y**2+z**2)
            part = simplify(part)
            result.append(part)
    return result


if __name__ == "__main__":
    def iter_m():
        counter = 0
        yield str(counter)
        while True:
            counter += 1
            yield "%i+" % counter
            yield "%i-" % counter

    x = Symbol("x")
    y = Symbol("y")
    z = Symbol("z")
    xyz = (x, y, z)
    labels = []
    print "solid_harmonics = ["
    for shell in xrange(5):
        im = iter_m()
        factor = sqrt(4*pi/(2*shell+1))
        for poly in get_solid_harmonics(shell, xyz):
            m = im.next()
            print "    lambda x,y,z: %s, # (%i,%s)" % (
                (poly*factor).evalf(20), shell, m
            )
            labels.append("(%i,%s)" % (shell, m))
    print "]"
    print "labels = %s" % labels

    print
    print "========  ", "="*150
    for shell in xrange(5):
        im = iter_m()
        factor = sqrt(4*pi/(2*shell+1))
        for poly in get_solid_harmonics(shell, xyz):
            m = im.next()
            print "%8s   :math:`%s`" % ("(%i,%s)" % (shell, m), latex(simplify(poly*factor)))
    print "========  ", "="*150
