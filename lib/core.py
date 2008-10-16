# HiPart is a tool to analyse molecular densities with the hirshfeld partitioning scheme
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

from molmod.data.periodic import periodic

import numpy, sys, os, shutil, pylab


__all__ = [
    "Error", "ProgressBar",
    "integrate_equi", "integrate", "cumul_integrate_equi", "cumul_integrate",
    "AtomFn", "AtomProfile", "AtomTable",
]


class Error(Exception):
    pass


class ProgressBar(object):
    def __init__(self, label, n, f=sys.stdout):
        self.i = 0
        self.n = n
        self.f = f
        self.f.write("%s: " % label)
        self.f.flush()

    def __call__(self):
        if self.i % 10 == 0 or self.i == self.n:
            self.f.write(" %i%% " % ((100*self.i)/self.n))
        else:
            self.f.write(".")
        if self.i == self.n:
            self.f.write("\n")
        self.f.flush()
        self.i += 1


noble_numbers = numpy.array([0,2,10,18,36,54,86,118])
core_sizes = dict((number, noble_numbers[noble_numbers<=number].max()) for number in periodic.yield_numbers())


def integrate_equi(u, g):
    s = CubicSpline(u, g)
    return s.integrate()

def integrate(r, f):
    u = numpy.log(r)
    g = f*r
    return integrate_equi(u,g)

def cumul_integrate_equi(u,g):
    s = CubicSpline(u, g)
    return s.cumul_integrate()

def cumul_integrate(r, f):
    u = numpy.log(r)
    g = f*r
    return cumul_integrate_equi(u,g)


class AtomFn(object):
    def __init__(self, rs, rhos, number, charge):
        self.number = number
        self.charge = charge
        self.num_elec = (number - charge)

        self.density = CubicSpline(rs, rhos)

        qs = -cumul_integrate(rs, rhos*rs**2*4*numpy.pi)
        #self.rhos_cumul = cumul_integrate(rs, rhos)
        #print number, self.rhos_cumul[-1], self.num_elec
        qs *= self.num_elec/qs[-1]
        self.charge = CubicSpline(rs, qs)

        vs = -cumul_integrate(1/rs[::-1], qs[::-1])[::-1]
        self._potential = CubicSpline(rs, vs)

    def potential(self, rs):
        vs = self._potential(rs)
        mask = rs>0#self._potential.x[-2]
        vs[mask] = -self.num_elec/rs[mask]
        return vs


class AtomProfile(object):
    def __init__(self, number, rs, records):
        self.number = number

        mask = reduce(
            (lambda x,y: x|y),
            (rhos > 0 for rhos in records.itervalues()),
            False
        )
        self.rs = rs[mask]
        #from molmod.units import angstrom
        #print number, self.rs[-1]/angstrom
        self.records = dict((charge, rhos[mask]) for charge, rhos in records.iteritems())

    def get_atom_fn(self, charge):
        if charge > self.number:
            raise Error("A negative number of electrons is not physical")

        max_charge = max(self.records)
        min_charge = min(self.records)
        if charge >= max_charge:
            if self.number > 1 and charge < 1:
                print "Warning: unsafe extrapolation (pos), number=%i, charge=%f" % (self.number, charge)
            num_elec = self.number - charge
            ref_elec = self.number - max_charge
            rhos = (self.records[max_charge]*num_elec)/ref_elec
        elif charge <= min_charge:
            print "Warning: unsafe extrapolation (neg), number=%i, charge=%f" % (self.number, charge)
            rhos = self.records[min_charge] + \
                (self.records[min_charge] - self.records[min_charge+1])*(min_charge-charge)
        else:
            high_charge = int(numpy.ceil(charge))
            low_charge = int(numpy.floor(charge))
            if low_charge == high_charge:
                low_charge -= 1
            low_ref = self.records[low_charge]
            high_ref = self.records[high_charge]
            rhos = high_ref + (low_ref - high_ref)*(high_charge - charge)

        result = AtomFn(self.rs, rhos, self.number, charge)
        #print "##Check:", -integrate(result.rs, 4*numpy.pi*result.rhos*result.rs**2)+self.number-charge, "##"
        return result

    def init_cusp_cutoff(self):
        if hasattr(self, "cusp_cutoff"):
            return
        rhos = self.records[0]
        qs = cumul_integrate(self.rs, rhos*self.rs**2*numpy.pi*4)
        i = qs.searchsorted([core_sizes[self.number]])[0]
        self.cusp_cutoff = self.rs[i]
        #print self.number, core_sizes[self.number], self.cusp_cutoff


class AtomTable(object):
    def __init__(self, filename):
        f = file(filename)
        line = f.next()
        self.rs = numpy.array([float(word) for word in line.split()[2:]])
        records = {}
        for line in f:
            words = line.split()
            number = int(words[1])
            charge = int(words[3])
            atom_rhos = numpy.array([float(word) for word in words[5:]])
            qmap = records.setdefault(number, {})
            qmap[charge] = atom_rhos
        f.close()
        self.records = {}
        for number, qmap in records.iteritems():
            self.records[number] = AtomProfile(number, self.rs, qmap)

    def init_cusp_cutoffs(self):
        for ad in self.records.itervalues():
            ad.init_cusp_cutoff()


