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


import numpy


__all__ = ["GaussianBasis"]


def get_shell_dof(shell_type):
    if shell_type >= 0:
        return 2*shell_type + 1
    elif shell_type == -1:
        return 8
    else:
        return (-shell_type+1)*(-shell_type+2)/2


class GaussianBasis(object):
    def __init__(self, molecule, shell_types, shell_map, num_primitives, ccoeffs, exponents):
        """
           Arguments:
            | ``molecule``  --  a Molecule object with atomic numbers and coordinates.
            | ``shell_types``  --  An array with shell types: 0 = S, 1 = P,
                                   2 = pure D, 3 = pure F, ..., -1 = SP,
                                   -2 = Cartesian D, -3 = Cartesian F, ...
            | ``shell_map``  --  An array with the atom index for each shell.
            | ``num_primitives``  --  The number of primitives in each shell.
            | ``ccoeffs``  --  The contraction coefficients of the
                                          primitives. For SP shells the number
                                          of contraction coefficients is twice
                                          the number of primitives, first n for
                                          S, second n for P.
            | ``exponents``  --  The exponents of the primitives.

           Standard convention for basis functions associated with a shell.

           The order of the pure shells is based on the order of real spherical
           harmonics: http://en.wikipedia.org/wiki/Table_of_spherical_harmonics
           First the +- linear combination of highest angular momentum, then
           the ++ combination of highest angular momentum, keep repeating and
           finally take angular momention zero (without making a linear
           combination). The order of the Cartesian shells is sorted
           alhpabetically. The SP shell type is S first, then P. Some examples:

           shell_type=0, S:
             0 -> 1
           shell_type=1, P:
             0 -> x
             1 -> y
             2 -> z
           shell_type=2, pure D:
             0 -> xy
             1 -> xx-yy
             2 -> xz
             3 -> yz
             4 -> zz
           shell_type=3, pure F:
             0 -> xxx-3xyy
             1 -> 3xxy-yyy
             2 -> xyz
             3 -> xxz-yyz
             4 -> xzz
             5 -> yzz
             6 -> zzz
           shell_type=-1, SP:
             0 -> 1
             1 -> x
             2 -> y
             3 -> z
           shell_type=-2, Cartesian D:
             0 -> xx
             1 -> xy
             2 -> xz
             3 -> yy
             4 -> yz
             5 -> zz
           shell_type=-3, Cartesian F:
             0 -> xxx
             1 -> xxy
             2 -> xxz
             3 -> xyy
             4 -> xyz
             5 -> xzz
             6 -> yyy
             7 -> yyz
             8 -> yzz
             9 -> zzz
        """
        self.molecule = molecule
        self.shell_types = shell_types
        self.shell_map = shell_map
        self.num_primitives = num_primitives
        self.ccoeffs = ccoeffs
        self.exponents = exponents
        # internal stuff
        self.num_dof = sum(get_shell_dof(shell_type) for shell_type in self.shell_types)
        self._prim_c_map = []
        self._prim_e_map = []
        counter_c = 0
        counter_e = 0
        for i in xrange(self.num_shells):
            if self.shell_types[i] == -1:
                size_c = 2*num_primitives[i]
                size_e = num_primitives[i]
            else:
                size_c = num_primitives[i]
                size_e = num_primitives[i]
            self._prim_c_map.append([counter_c, size_c])
            self._prim_e_map.append([counter_e, size_e])
            counter_c += size_c
            counter_e += size_e

    num_shells = property(lambda self: len(self.shell_types))

    @classmethod
    def from_fchk(cls, fchk):
        shell_types = fchk.fields["Shell types"]
        shell_map = fchk.fields["Shell to atom map"] - 1
        num_primitives = fchk.fields["Number of primitives per shell"]
        ccoeffs_level1 = fchk.fields["Contraction coefficients"]
        ccoeffs_level2 = fchk.fields.get("P(S=P) Contraction coefficients")
        exponents = fchk.fields["Primitive exponents"]

        ccoeffs = []
        counter = 0
        for i, n in enumerate(num_primitives):
            ccoeffs.append(ccoeffs_level1[counter:counter+n])
            if shell_types[i] == -1:
                ccoeffs.append(ccoeffs_level2[counter:counter+n])
        ccoeffs = numpy.concatenate(ccoeffs)

        return cls(fchk.molecule, shell_types, shell_map, num_primitives, ccoeffs, exponents)

    def call_gint1(self, gint1_fn, weights, point):
        output = numpy.zeros(1, float)
        retcode = gint1_fn(
            weights, output, point, self.molecule.coordinates, self.shell_types,
            self.shell_map, self.num_primitives, self.ccoeffs, self.exponents
        )
        if retcode != 0:
            raise RunTimeError("Something went wrong when calling %s. Got retcode=%i." % (gint1_fn, retcode))
        return output[0]
