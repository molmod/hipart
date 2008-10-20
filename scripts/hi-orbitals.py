#!/usr/bin/python
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


from hipart.context import Context
from hipart.atoms import AtomTable
from hipart.opts import add_hirshi_options

from optparse import OptionParser
import os, numpy


# parse arguments ...
parser = OptionParser("""%prog [options] atoms.txt gaussian.fchk

%prog partitions the density matrices with iterative hirshfeld weights.
""")
add_hirshi_options(parser)
(options, args) = parser.parse_args()
if len(args) != 2:
    parser.error("Expecting two arguments: atoms.txt gaussian.fchk")
atom_fn, fchk_fn = args

# Do the work. Where possible, the intermediate results from scripts that ran
# previously, are recycled.
context = Context(AtomTable(atom_fn), fchk_fn, options)
context.cache.do_atom_matrices()
context.clean()

