# -*- coding: utf-8 -*-
# HiPart is a program to analyze the electronic structure of molecules with
# fuzzy-atom partitioning methods.
# Copyright (C) 2007 - 2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
#--


import os, numpy, shutil


__all__ = ["Work"]


class Work(object):
    def __init__(self, directory=None, do_clean=False):
        self.directory = directory
        self.do_clean = do_clean
        self.clean()
        self.create()

    active = property(lambda self: self.directory is not None)

    def create(self):
        if self.active and not os.path.isdir(self.directory):
            os.makedirs(self.directory)

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
            array.tofile(filename)

    def clean(self):
        if self.do_clean and self.active and os.path.isdir(self.directory):
            shutil.rmtree(self.directory)
