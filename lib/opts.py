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


from hipart.lebedev_laikov import grid_fns
from hipart.cache import cache_classes, ParseError
from hipart.context import Context

from optparse import OptionParser


__all__ = ["parse_command_line"]

usage_template = """%%prog [options] gaussian.fchk scheme [scheme parameters]

%s
"""

def parse_command_line(script_usage):
    parser = OptionParser(usage_template % script_usage)
    parser.add_option(
        "--density",
        help="The density field to use from the gaussian fchk file (scf, mp2, mp3, "
        "...). If not given, the program will guess it based on the level of "
        "theory."
    )
    parser.add_option(
        "-l", "--lebedev", default=110, type='int',
        help="The number of grid points for the atomic grids. "
        "[default=%default]. Select from: " + (", ".join(str(i) for i in sorted(grid_fns)))
    )
    parser.add_option(
        "-m", "--mol-lebedev", default=50, type='int',
        help="The number of grid points for the molecular grids. "
        "[default=%default]. See --lebedev for the supported grid sizes."
    )
    parser.add_option(
        "-r", "--reference", default=None,
        help="Refer to the checkpoint af a reference state. The output will "
        "relative changes instead absolute values."
    )
    parser.add_option(
        "-c", "--clean", default=1, type='int',
        help="Degree of cleaning in the workdir after the computations are done. "
        "Files in the workdir can be reused by other scripts, which reduces the "
        "computational cost. This should be a number from 0 to 4. "
        "[default=%default] 0: No cleaning. 1: Remove text files. 2: Also remove "
        "binary files. 3: Remove the entire workdir."
    )
    parser.add_option(
        "-n", "--no-fix-total-charge", dest="fix_total_charge", default="True",
        action="store_false", help="Do not correct the total charge."
    )
    parser.add_option(
        "-t", "--threshold", default=1e-4, type='float',
        help="When the maximum change in the charges drops below this threshold "
        "value, the iteration stops. [default=%default]"
    )
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error("Expecting at least two arguments.")
    fchk_fn, scheme_name = args[:2]

    context = Context(fchk_fn, options)

    CacheClass = cache_classes.get(scheme_name)
    if CacheClass is None:
        parser.error("The scheme must be one of: %s" % (" ".join(sorted(cache_classes))))
    try:
        cache = CacheClass.new_from_args(context, args[2:])
    except ParseError, e:
        parser.error(str(e))

    return context, cache


