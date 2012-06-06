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



import os, tempfile, shutil


__all__ = ["setup_gaussian"]


def setup_gaussian(fchk_name):
    tmpdir = tempfile.mkdtemp("hipart")
    if not os.path.isdir("input"):
        raise IOError("Input directory with test files is not present")
    fn_fchk = os.path.join(tmpdir, "gaussian.fchk")
    shutil.copy("input/%s.fchk" % fchk_name, fn_fchk)
    return tmpdir, fn_fchk
