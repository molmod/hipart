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






from hipart.lebedev_laikov import grid_fns
from hipart.schemes import scheme_classes, ParseError
from hipart.context import Context, ContextError

from optparse import OptionParser, OptionGroup
import sys


__all__ = ["parse_command_line"]

usage_template = """%%prog [options] gaussian.fchk scheme [scheme parameters]

%s

The file gaussian.fchk is a formatted checkpoint file from a Gaussian
computation. To obtain this file, add the following line on top of a Gaussian
com-file (before running the job)

%%chk=gaussian.chk

After the Gaussian computation transform this binary checkpoint file into
a text file with the ``formchk`` program of the Gaussian software suite:

formchk gaussian.chk gaussian.fchk

Partitioning schemes:

%s
"""

def parse_command_line(script_usage):
    scheme_usage = "\n".join(scheme.usage for name, scheme in sorted(scheme_classes.iteritems()))
    parser = OptionParser(usage_template % (script_usage, scheme_usage))
    parser.add_option(
        "-l", "--lebedev", default=110, type='int',
        help="The number of grid points for the atomic grids. "
        "[default=%default]. Select from: " + (", ".join(str(i) for i in sorted(grid_fns)))
    )
    parser.add_option(
        "-c", "--clean", default=False, action='store_true', dest='do_clean',
        help="Remove the workdir before and after the computation."
    )
    parser.add_option(
        "--no-work", default=True, action='store_false', dest='do_work',
        help="Do not save intermediate results in work directory for later reuse."
    )
    parser.add_option(
        "--no-output", default=True, action='store_false', dest='do_output',
        help="Do not write any output to text files."
    )
    parser.add_option(
        "-q", "--quiet", default=True, action='store_false', dest='verbose',
        help="Do not write any screen output."
    )
    parser.add_option(
        "--no-fix-total-charge", dest="fix_total_charge", default=True,
        action="store_false", help="Do not correct the total charge."
    )
    parser.add_option(
        "--no-random", default=True, action='store_false', dest='do_random',
        help="Do not randomly rotate angular grids."
    )
    parser.add_option(
        "--save-mem", action="store_true", default=False,
        help="Try to be less memory hungry at the expense of a little efficiency."
    )
    group = OptionGroup(parser, "Specific options for the iterative partitioning schemes")
    group.add_option(
        "-t", "--threshold", default=1e-4, type='float',
        help="When the maximum change in the charges drops below this threshold "
        "value, the iteration stops. [default=%default]"
    )
    group.add_option(
        "--max-iter", default=500, type='int',
        help="Maximum number of iterations in self-consistent procedures. "
        "[default=%default]"
    )
    parser.add_option_group(group)

    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error("Expecting at least two arguments.")
    fchk_fn, scheme_name = args[:2]

    context = Context(fchk_fn, options)

    SchemeClass = scheme_classes.get(scheme_name)
    if SchemeClass is None:
        parser.error("The scheme must be one of: %s" % (" ".join(sorted(scheme_classes))))
    try:
        scheme = SchemeClass.new_from_args(context, args[2:])
    except ParseError, e:
        print str(e)
        sys.exit(-1)
    except ContextError, e:
        print str(e)
        sys.exit(-1)

    return context, scheme
