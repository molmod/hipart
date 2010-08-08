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


import os, numpy


__all__ = ["Work"]


class Work(object):
    def __init__(self, directory=None, do_clean=True):
        self.directory = directory
        self.do_clean = do_clean
        if self.active and not os.path.isdir(self.directory):
            os.makedirs(self.directory)

    active = property(lambda self: self.directory is not None)

    def load(self, name, shape=None):
        if self.active:
            filename = os.path.join(self.directory, name + ".bin")
            if os.path.isfile(filename):
                array = numpy.fromfile(filename)
                if shape is not None:
                    array = array.reshape(shape)
                return array

    def dump(self, name, array, ignore=False):
        if self.active:
            filename = os.path.join(self.directory, name + ".bin")
            if os.path.isfile(filename):
                if not ignore:
                    raise ValueError("The binary file is already present in the work directory.")
            else:
                array.tofile(filename)

    def clean(self):
        if self.do_clean and self.active and os.path.isdir(self.workdir):
            shutil.rmtree(self.workdir)
