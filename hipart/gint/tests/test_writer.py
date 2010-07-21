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


from hipart.gint.writer import *


def test_cartesian_powers():
    assert(list(get_cartesian_powers(0))==[(0,0,0)])
    assert(list(get_cartesian_powers(1))==[(1,0,0), (0,1,0), (0,0,1)])
    assert(list(get_cartesian_powers(2))==[(2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2)])

def test_norms():
    alpha = Symbol("alpha")
    a = simplify(mypowsimp(get_cartesian_wfn_norm(alpha, 0, 0, 0)))
    b = simplify(mypowsimp(get_pure_wfn_norm(alpha, 0)))
    print a
    print b
    assert(a==b)
    a = simplify(mypowsimp(get_cartesian_wfn_norm(alpha, 1, 0, 0)))
    b = simplify(mypowsimp(get_pure_wfn_norm(alpha, 1)))
    print a
    print b
    print b/a
    assert(a==b)
