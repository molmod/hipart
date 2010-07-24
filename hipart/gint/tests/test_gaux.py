#!/usr/bin/env python
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

from hipart.gint.gaux import gaux

import numpy
from scipy.special import erf


def test_gaux():
    for t in 0.1, 0.5, 0.7, 1.0, 2.0, 4.0:
        u1 = gaux(t, 0)
        u2 = numpy.sqrt(numpy.pi/t)/2*erf(numpy.sqrt(t))
