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

from hipart.gint.tests.utils import setup_fchk, h_sto3g_fchk, hf_sto3g_fchk, \
    o2_cc_pvtz_pure_fchk, o2_cc_pvtz_cart_fchk

from hipart.gint.basis import GaussianBasis, get_shell_dof

from molmod.io import FCHKFile

import shutil, numpy


def test_hf_sto3g_num():
    tmpdir, fn_fchk = setup_fchk(hf_sto3g_fchk)
    fchk = FCHKFile(fn_fchk)
    basis = GaussianBasis.from_fchk(fchk)
    assert(basis.num_shells==3)
    assert(basis.num_dof==6)


def test_h_sto3g_num():
    tmpdir, fn_fchk = setup_fchk(h_sto3g_fchk)
    fchk = FCHKFile(fn_fchk)
    basis = GaussianBasis.from_fchk(fchk)
    assert(basis.num_shells==1)
    assert(basis.num_dof==1)


def test_o2_cc_pvtz_pure_num():
    tmpdir, fn_fchk = setup_fchk(o2_cc_pvtz_pure_fchk)
    fchk = FCHKFile(fn_fchk)
    basis = GaussianBasis.from_fchk(fchk)
    assert(basis.num_shells==20)
    print basis.num_dof
    assert(basis.num_dof==60)


def test_o2_cc_pvtz_cart_num():
    tmpdir, fn_fchk = setup_fchk(o2_cc_pvtz_cart_fchk)
    fchk = FCHKFile(fn_fchk)
    basis = GaussianBasis.from_fchk(fchk)
    assert(basis.num_shells==20)
    print basis.num_dof
    assert(basis.num_dof==70)


def test_shell_dof():
    assert(get_shell_dof(-3)==7)
    assert(get_shell_dof(-2)==5)
    assert(get_shell_dof(-1)==4)
    assert(get_shell_dof( 0)==1)
    assert(get_shell_dof( 1)==3)
    assert(get_shell_dof( 2)==6)
    assert(get_shell_dof( 3)==10)
