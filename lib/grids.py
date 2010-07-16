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


from tools import relsymlink
import numpy, os


__all__ = ["Grid"]


class Grid(object):
    def __init__(self, prefix, points, dump=True):
        self.prefix = prefix
        self.points = points
        if dump:
            fn = "%s.bin" % prefix
            if os.path.isfile(fn):
                raise ValueError("The binary file is already present in the work directory.")
            else:
                points.tofile(fn)

    @classmethod
    def from_prefix(cls, prefix):
        fn = "%s.bin" % prefix
        if os.path.isfile(fn):
            points = numpy.fromfile(fn).reshape((-1,3))
            return cls(prefix, points, False)
        else:
            return None

    def get_fn(self, suffix, ext="bin"):
        return "%s_%s.%s" % (self.prefix, suffix, ext)

    def load(self, suffix):
        fn = self.get_fn(suffix)
        if os.path.isfile(fn):
            return numpy.fromfile(fn)
        else:
            return None

    def dump(self, suffix, array):
        fn = self.get_fn(suffix)
        if os.path.isfile(fn):
            raise ValueError("The binary file is already present in the work directory.")
        else:
            return array.tofile(fn)


