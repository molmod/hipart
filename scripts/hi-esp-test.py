#!/usr/bin/python
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






from hipart.opts import parse_command_line

usage = """%prog test charges and dipoles on ESP grid

These atomic charges and dipoles are computed like in the hi-charges.py and
hi-dipoles.py program. Then a ESP costfunction is constructed and it is tested
how well the partitioned charges and dipoles reproduce the ESP."""

def add_extra_options(group):
    group.add_option(
        "-m", "--mol-lebedev", default=50, type='int',
        help="The number of grid points for the molecular grids. "
        "[default=%default]. See --lebedev for the supported grid sizes."
    )

context, cache = parse_command_line(usage, add_extra_options)
cache.do_esp_test()
context.clean()
