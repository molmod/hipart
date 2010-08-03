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


from hipart.log import log
from hipart.spline import CubicSpline

import numpy, os


__all__ = [
    "AtomProfile", "AtomTable",
]


class AtomProfile(object):
    def __init__(self, number, rs, records):
        self.number = number

        mask = reduce(
            (lambda x,y: x|y),
            (rhos > 0 for rhos in records.itervalues()),
            False
        )
        self.rs = rs[mask]
        self.records = dict((charge, rhos[mask]) for charge, rhos in records.iteritems())

    def get_atom_fn(self, charge):
        if charge > self.number:
            raise ValueError("A negative number of electrons is not physical")

        max_charge = max(self.records)
        min_charge = min(self.records)
        if charge >= max_charge:
            if self.number > 1 and charge < 1:
                log("Warning: unsafe extrapolation (pos), number=%i, charge=%f" % (self.number, charge))
            num_elec = self.number - charge
            ref_elec = self.number - max_charge
            rhos = (self.records[max_charge]*num_elec)/ref_elec
        elif charge <= min_charge:
            log("Warning: unsafe extrapolation (neg), number=%i, charge=%f" % (self.number, charge))
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

        result = CubicSpline(self.rs, rhos)
        #print "##Check:", -integrate(result.rs, 4*numpy.pi*result.rhos*result.rs**2)+self.number-charge, "##"
        return result


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
