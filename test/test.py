#!/usr/bin/env python
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


import sys, os, glob, shutil
if os.path.isdir("../build"):
    shutil.rmtree("../build")
os.system("cd ../; python setup.py build")
sys.path.insert(0, "../build/lib")
os.system("cd ../ext/; python setup.py build")
for filename in glob.glob("../ext/build/lib*/hipart/*"):
    shutil.copy(filename, "../build/lib/hipart")

if not os.path.isdir("output"):
    os.mkdir("output")

import unittest
from integrate import *
from grid import *
from spline import *
from tools import *
unittest.main()




