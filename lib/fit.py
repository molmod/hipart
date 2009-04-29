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


import numpy


__all__ = ["ChargeDipoleCostFunction", "ESPCostFunction"]


class ChargeDipoleCostFunction(object):
    def __init__(self, N, design_matrix, expected_values, weights, total_charge):
        self.N = N # number of atoms
        self.design_matrix = design_matrix*weights.reshape((-1,1))
        self.expected_values = expected_values*weights
        self.weights = weights
        self.weights_sqsum = (weights**2).sum()

        self.A = numpy.dot(self.design_matrix.transpose(), self.design_matrix)
        self.B = numpy.dot(self.design_matrix.transpose(), self.expected_values)
        self.C = numpy.dot(self.expected_values.transpose(), self.expected_values)
        self.rms = numpy.sqrt(self.C/self.weights_sqsum)

        Ap = numpy.zeros((len(self.A)+1,len(self.A)+1), float)
        Ap[:-1,:-1] = self.A
        Ap[:-1,-1] = 1
        Ap[-1,:-1] = 1
        Bp = numpy.zeros(len(self.B)+1, float)
        Bp[:-1] = self.B
        Bp[-1] = total_charge
        self.solution = numpy.linalg.solve(Ap, Bp)[:-1]

        A_evals = numpy.linalg.eigvalsh(self.A)
        self.condition_number = A_evals[-1]/A_evals[0]

    def full_vector(self, charges, dipoles=None):
        result = numpy.zeros(4*self.N, float)
        if charges is not None:
            result[:self.N] = charges
        if dipoles is not None:
            result[self.N:] = dipoles.ravel()
        return result

    def rmsd(self, charges, dipoles=None):
        full = self.full_vector(charges, dipoles)
        return numpy.sqrt((numpy.dot(full, numpy.dot(self.A, full) - 2*self.B) + self.C)/self.weights_sqsum)

    def model_rms(self, charges, dipoles=None):
        full = self.full_vector(charges, dipoles)
        return numpy.sqrt(numpy.dot(full, numpy.dot(self.A, full))/self.weights_sqsum)

    def correlation(self, charges, dipoles=None):
        # pearson's r
        full = self.full_vector(charges, dipoles)
        return numpy.dot(full, self.B)/numpy.sqrt(self.C*numpy.dot(full, numpy.dot(self.A, full)))

    def write_to_file(self, filename):
        f = file(filename, "w")
        print >> f, "number of atoms:", self.N
        print >> f, "A: Quadratic part of the ESP Cost function"
        for row in self.A:
            print >> f, " ".join("% 15.10e" % value for value in row)

        print >> f, "B: Linear part of the ESP Cost function"
        print >> f, " ".join("% 15.10e" % value for value in self.B)

        print >> f, "C: Constant term of the ESP Cost function"
        print >> f, "% 15.10e" % self.C

        print >> f, "W: Sum of the squared weights"
        print >> f, "% 15.10e" % self.weights_sqsum
        f.close()


class ESPCostFunction(ChargeDipoleCostFunction):
    def __init__(self, coordinates, grid_points, grid_weights, densities, potentials, total_charge):
        self.grid_points = grid_points
        self.grid_weights = grid_weights
        self.densities = densities
        self.potentials = potentials

        N = len(coordinates)
        expected_values = potentials.copy()
        design_matrix = numpy.zeros((len(expected_values), 4*N), float)
        self.distances = {}
        for i in xrange(N):
            deltas = grid_points - coordinates[i]
            distances = numpy.sqrt((deltas**2).sum(axis=1))
            self.distances[i] = distances
            design_matrix[:,i] = 1/distances
            design_matrix[:,N+3*i:N+3*i+3] = deltas/distances.reshape((-1,1))**3
        weights = grid_weights * numpy.exp(-densities/1e-5)
        ChargeDipoleCostFunction.__init__(self, N, design_matrix, expected_values, weights, total_charge)


