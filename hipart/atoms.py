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
from hipart.spline import CubicSpline, get_rgrid_from_description

import numpy, os


__all__ = [
    "AtomProfile", "AtomTable",
]


class AtomProfile(object):
    def __init__(self, number, rgrid, records):
        self.number = number
        self.rgrid = rgrid
        self.records = records
        self.max_population = max(self.records)
        self.min_population = min(self.records)

    def get_atom_fn(self, population=None):
        if population is None:
            population = float(self.number)
        else:
            population = float(population)
        if population < 0:
            raise ValueError("A negative number of electrons is not physical")

        if population < self.min_population:
            if self.number > 1 and population < self.min_population:
                log("Warning: unsafe extrapolation (below), number=%i, population=%f" % (self.number, population))
            ratio = population/self.min_population
            rhos = self.records[self.min_population]*ratio
        elif population > self.max_population:
            log("Warning: unsafe extrapolation (above), number=%i, population=%f" % (self.number, population))
            ratio = population/self.max_population
            rhos = self.records[self.max_population]*ratio
        else:
            high_population = int(numpy.ceil(population))
            low_population = int(numpy.floor(population))
            if low_population == high_population:
                rhos = self.records[low_population]
            else:
                low_ref = self.records[low_population]
                high_ref = self.records[high_population]
                if len(high_ref) > len(low_ref):
                    rhos = high_ref*(population - low_population)
                    rhos[:len(low_ref)] += low_ref*(high_population - population)
                else:
                    rhos = low_ref*(high_population - population)
                    rhos[:len(high_ref)] += high_ref*(population - low_population)

        result = CubicSpline(self.rgrid.rs[:len(rhos)], rhos)
        return result


class AtomTable(object):
    def __init__(self, filename):
        f = file(filename)
        line = f.next()
        self.rgrid = get_rgrid_from_description(line.strip())
        records = {}
        for line in f:
            words = line.split()
            number = int(words[0])
            charge = int(words[1])
            population = number - charge
            atom_rhos = numpy.array([float(word) for word in words[2:]])
            nmap = records.setdefault(number, {})
            nmap[population] = atom_rhos
        f.close()
        self.records = {}
        for number, nmap in records.iteritems():
            self.records[number] = AtomProfile(number, self.rgrid, nmap)
