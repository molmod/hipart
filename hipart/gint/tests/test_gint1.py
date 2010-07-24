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
    o2_cc_pvtz_cart_fchk, o2_cc_pvtz_pure_fchk

from hipart.gint.basis import GaussianBasis
from hipart.gint.gint1_fn import gint1_fn_basis, gint1_fn_dmat
from hipart.gint.ctools import reorder_density_matrix

from molmod.io import FCHKFile
from molmod import angstrom

import numpy, shutil


ref_data_h_sto3g_orb0 = numpy.array([
    [0.0, 0.0, 0.0, 0.628246876936],
    [0.0, 0.0, 0.1, 0.589959208052],
    [0.0, 0.0, 0.2, 0.496266720258],
    [0.0, 0.0, 0.3, 0.39034158533],
    [0.0, 0.0, 0.4, 0.302441507174],
    [0.0, 0.0, 0.5, 0.238254435193],
    [0.0, 0.0, 0.6, 0.190669678643],
    [0.0, 0.0, 0.7, 0.152706130579],
    [0.0, 0.0, 0.8, 0.121214786138],
    [0.0, 0.0, 0.9, 0.095298026935],
    [0.0, 0.0, 1.0, 0.074532634931],
    [0.0, 0.0, 1.1, 0.058315779579],
    [0.0, 0.0, 1.2, 0.045857230231],
    [0.0, 0.0, 1.3, 0.036331335546],
    [0.0, 0.0, 1.4, 0.029000669294],
    [0.0, 0.0, 1.5, 0.023275698483],
    [0.0, 0.0, 1.6, 0.018722540918],
    [0.0, 0.0, 1.7, 0.015040532334],
    [0.0, 0.0, 1.8, 0.012028520273],
    [0.0, 0.0, 1.9, 0.009552219886],
])


def test_orb0_h_sto3g():
    tmpdir, fn_fchk = setup_fchk(h_sto3g_fchk)
    fchk = FCHKFile(fn_fchk)
    shutil.rmtree(tmpdir)
    basis = GaussianBasis.from_fchk(fchk)
    weights = fchk.fields["Alpha MO coefficients"][:basis.num_dof]
    points = ref_data_h_sto3g_orb0[:,:3]
    fns = basis.call_gint1(gint1_fn_basis, weights, points*angstrom)
    assert(abs(fns-ref_data_h_sto3g_orb0[:,3]).max() < 1e-10)


ref_data_hf_sto3g_orb0 = numpy.array([
    [0.0, 0.0, 0.0, 2.789278788847],
    [0.0, 0.0, 0.1, 11.520724973181],
    [0.0, 0.0, 0.2, 2.857261533838],
    [0.0, 0.0, 0.3, 0.552008457019],
    [0.0, 0.0, 0.4, 0.118807434542],
    [0.0, 0.0, 0.5, 0.020878797978],
    [0.0, 0.0, 0.6, 0.006053066021],
    [0.0, 0.0, 0.7, 0.003703387901],
    [0.0, 0.0, 0.8, 0.002632278844],
    [0.0, 0.0, 0.9, 0.001876910771],
    [0.0, 0.0, 1.0, 0.001322579924],
    [0.0, 0.0, 1.1, 0.000913262609],
    [0.0, 0.0, 1.2, 0.000613437551],
    [0.0, 0.0, 1.3, 0.000398715852],
    [0.0, 0.0, 1.4, 0.000249908777],
    [0.0, 0.0, 1.5, 0.000150676115],
    [0.0, 0.0, 1.6, 0.000087174479],
    [0.0, 0.0, 1.7, 0.000048237147],
    [0.0, 0.0, 1.8, 0.000025393271],
    [0.0, 0.0, 1.9, 0.000012599175],
])

def test_orb0_hf_sto3g():
    tmpdir, fn_fchk = setup_fchk(hf_sto3g_fchk)
    fchk = FCHKFile(fn_fchk)
    shutil.rmtree(tmpdir)
    basis = GaussianBasis.from_fchk(fchk)
    weights = fchk.fields["Alpha MO coefficients"][:basis.num_dof]
    points = ref_data_hf_sto3g_orb0[:,:3]
    fns = basis.call_gint1(gint1_fn_basis, weights, points*angstrom)
    assert(abs(fns-ref_data_hf_sto3g_orb0[:,3]).max() < 1e-10)


ref_data_o2_cc_pvtz_cart_orb0 = numpy.array([
    [0.0, 0.0, 0.0, 0.007491177473],
    [0.0, 0.0, 0.2, 0.043202323599],
    [0.0, 0.0, 0.4, 0.657012617121],
    [0.0, 0.0, 0.6, 6.131445881796],
    [0.0, 0.0, 0.8, 0.345132257671],
    [0.0, 0.0, 1.0, 0.022026771467],
    [0.0, 0.0, 1.2, 0.002082949266],
    [0.0, 0.0, 1.4, 0.000179713422],
    [0.0, 0.0, 1.6, 0.000063698692],
    [0.0, 0.0, 1.8, 0.000048544679],
    [-0.560286181718, 0.909915321675, 0.41252418403, 0.000086938833],
    [1.24430316424, -1.4500487541, 0.099047495872, -0.000016512327],
    [0.169854001784, -0.44672091063, -0.510747144295, 0.010301854004],
    [-1.90578665178, 0.24668484101, 1.04268268033, -0.000012757029],
    [-0.628611591276, 0.806681447401, -0.381539890553, 0.000091028151],
    [-0.670647186735, 1.19432378315, 0.719831575024, 0.000023666149],
    [1.53885489642, 2.06941939659, 0.811533561058, -0.000001442113],
    [0.523958857852, -0.305479070895, 0.337511300396, 0.001634861602],
    [0.266993506916, 0.470266561219, 0.655807944137, 0.004748956154],
    [0.448869071646, 0.94063466163, 1.98016240483, -0.000022014793],
])

def test_orb0_o2_cc_pvtz_cart():
    tmpdir, fn_fchk = setup_fchk(o2_cc_pvtz_cart_fchk)
    fchk = FCHKFile(fn_fchk)
    shutil.rmtree(tmpdir)
    basis = GaussianBasis.from_fchk(fchk)
    weights = fchk.fields["Alpha MO coefficients"][:basis.num_dof]
    weights = weights[basis.g03_permutation]
    points = ref_data_o2_cc_pvtz_cart_orb0[:,:3]
    fns = basis.call_gint1(gint1_fn_basis, weights, points*angstrom)
    assert(abs(fns-ref_data_o2_cc_pvtz_cart_orb0[:,3]).max() < 1e-10)


ref_data_o2_cc_pvtz_pure_orb0 = numpy.array([
    [0.0, 0.0, 0.0, 0.007424022071],
    [0.0, 0.0, 0.2, 0.043139399479],
    [0.0, 0.0, 0.4, 0.657048305371],
    [0.0, 0.0, 0.6, 6.132762906955],
    [0.0, 0.0, 0.8, 0.345173467373],
    [0.0, 0.0, 1.0, 0.021993579251],
    [0.0, 0.0, 1.2, 0.002220896203],
    [0.0, 0.0, 1.4, 0.000188868498],
    [0.0, 0.0, 1.6, -0.000036771101],
    [0.0, 0.0, 1.8, -0.000002706907],
    [-0.560286181718, 0.909915321675, 0.41252418403, 0.000043390084],
    [1.24430316424, -1.4500487541, 0.099047495872, 0.000025605707],
    [0.169854001784, -0.44672091063, -0.510747144295, 0.010272778148],
    [-1.90578665178, 0.24668484101, 1.04268268033, 0.000012360357],
    [-0.628611591276, 0.806681447401, -0.381539890553, 0.00003585536],
    [-0.670647186735, 1.19432378315, 0.719831575024, 0.000060676876],
    [1.53885489642, 2.06941939659, 0.811533561058, 0.000001339715],
    [0.523958857852, -0.305479070895, 0.337511300396, 0.001675874092],
    [0.266993506916, 0.470266561219, 0.655807944137, 0.004772119841],
    [0.448869071646, 0.94063466163, 1.98016240483, 0.000013101125],
])

def test_orb0_o2_cc_pvtz_pure():
    tmpdir, fn_fchk = setup_fchk(o2_cc_pvtz_pure_fchk)
    fchk = FCHKFile(fn_fchk)
    shutil.rmtree(tmpdir)
    basis = GaussianBasis.from_fchk(fchk)
    weights = fchk.fields["Alpha MO coefficients"][:basis.num_dof]
    weights = weights[basis.g03_permutation]
    points = ref_data_o2_cc_pvtz_pure_orb0[:,:3]
    fns = basis.call_gint1(gint1_fn_basis, weights, points*angstrom)
    assert(abs(fns-ref_data_o2_cc_pvtz_pure_orb0[:,3]).max() < 1e-10)


def test_dens_h_sto3g():
    tmpdir, fn_fchk = setup_fchk(h_sto3g_fchk)
    fchk = FCHKFile(fn_fchk)
    shutil.rmtree(tmpdir)
    basis = GaussianBasis.from_fchk(fchk)
    weights = fchk.fields["Alpha MO coefficients"][:basis.num_dof]
    points = ref_data_h_sto3g_orb0[:,:3]
    orb0 = basis.call_gint1(gint1_fn_basis, weights, points*angstrom)
    num_dmat = (basis.num_dof*(basis.num_dof+1))/2
    dmat = fchk.fields["Total SCF Density"][:num_dmat]
    density = basis.call_gint1(gint1_fn_dmat, dmat, points*angstrom)
    assert(abs(orb0**2-density).max() < 1e-10)


def test_dens_hf_sto3g():
    tmpdir, fn_fchk = setup_fchk(hf_sto3g_fchk)
    fchk = FCHKFile(fn_fchk)
    shutil.rmtree(tmpdir)
    basis = GaussianBasis.from_fchk(fchk)
    num_dmat = (basis.num_dof*(basis.num_dof+1))/2

    points = ref_data_hf_sto3g_orb0[:,:3]

    dmat = fchk.fields["Total SCF Density"][:num_dmat]
    reorder_density_matrix(dmat, basis.g03_permutation)
    density = basis.call_gint1(gint1_fn_dmat, dmat, points*angstrom)

    num_alpha = fchk.fields["Number of alpha electrons"]
    expected_density = 0.0
    start = 0
    for i in xrange(num_alpha):
        end = start + basis.num_dof
        weights = fchk.fields["Alpha MO coefficients"][start:end]
        start = end
        orb = basis.call_gint1(gint1_fn_basis, weights, points*angstrom)
        expected_density += 2*(orb**2)

    assert(abs(density-expected_density).max() < 1e-6)


def test_dens_o2_cc_pvtz_cart():
    tmpdir, fn_fchk = setup_fchk(o2_cc_pvtz_cart_fchk)
    fchk = FCHKFile(fn_fchk)
    shutil.rmtree(tmpdir)
    basis = GaussianBasis.from_fchk(fchk)
    num_dmat = (basis.num_dof*(basis.num_dof+1))/2

    points = ref_data_o2_cc_pvtz_cart_orb0[:,:3]

    dmat = fchk.fields["Total SCF Density"][:num_dmat]
    reorder_density_matrix(dmat, basis.g03_permutation)
    density = basis.call_gint1(gint1_fn_dmat, dmat, points*angstrom)

    num_alpha = fchk.fields["Number of alpha electrons"]
    expected_density = 0.0
    start = 0
    for i in xrange(num_alpha):
        end = start + basis.num_dof
        weights = fchk.fields["Alpha MO coefficients"][start:end]
        start = end
        orb = basis.call_gint1(gint1_fn_basis, weights, points*angstrom)
        expected_density += 2*(orb**2)

    assert(abs(density-expected_density).max() < 1e-6)


def test_dens_o2_cc_pvtz_pure():
    tmpdir, fn_fchk = setup_fchk(o2_cc_pvtz_pure_fchk)
    fchk = FCHKFile(fn_fchk)
    shutil.rmtree(tmpdir)
    basis = GaussianBasis.from_fchk(fchk)
    num_dmat = (basis.num_dof*(basis.num_dof+1))/2

    points = ref_data_o2_cc_pvtz_pure_orb0[:,:3]

    dmat = fchk.fields["Total SCF Density"][:num_dmat]
    reorder_density_matrix(dmat, basis.g03_permutation)
    density = basis.call_gint1(gint1_fn_dmat, dmat, points*angstrom)

    num_alpha = fchk.fields["Number of alpha electrons"]
    expected_density = 0.0
    start = 0
    for i in xrange(num_alpha):
        end = start + basis.num_dof
        weights = fchk.fields["Alpha MO coefficients"][start:end]
        start = end
        orb = basis.call_gint1(gint1_fn_basis, weights, points*angstrom)
        expected_density += 2*(orb**2)

    assert(abs(density-expected_density).max() < 1e-6)
