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


from hipart.tests.utils import iter_hf_sto3g_gaussian_schemes, \
    iter_oh1_sto3g_gaussian_schemes, iter_oh2_sto3g_gaussian_schemes, \
    iter_h_sto3g_gaussian_schemes
from hipart.gint import dmat_to_full

import numpy, os
from nose.plugins.skip import SkipTest


def test_compute_atgrid_atweights():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_compute_atgrid_atweights, scheme

def check_compute_atgrid_atweights(scheme):
    scheme.do_atgrids()
    scheme._prepare_atweights()
    h0 = scheme._compute_atweights(scheme.atgrids[0], 0)
    h1 = scheme._compute_atweights(scheme.atgrids[0], 1)
    error = abs(h1+h0-1).max()
    assert(error < 1e-10)


def test_hf_charges():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_hf_charges, scheme

def check_hf_charges(scheme):
    expected = {
        "hirsh": numpy.array([-0.13775422, 0.13775422]),
        "hirshi": numpy.array([-0.19453209, 0.19453209]),
        "isa": numpy.array([-0.20842675, 0.20842675]),
        "becke": numpy.array([-0.20037239, 0.20037239]),
    }
    scheme.do_charges()
    assert(abs(scheme.charges - expected[scheme.prefix]).max() < 1e-4)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_charges.txt" % scheme.prefix)))


def test_oh1_charges():
    for scheme in iter_oh1_sto3g_gaussian_schemes():
        yield check_oh1_charges, scheme

def check_oh1_charges(scheme):
    expected = {
        "hirsh": numpy.array([-0.11261, 0.11261]),
        "hirshi": numpy.array([-0.17905023, 0.17905023]),
        "isa": numpy.array([-0.19495561, 0.19495561]),
    }
    scheme.do_charges()
    assert(abs(scheme.charges - expected[scheme.prefix]).max() < 1e-3)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_charges.txt" % scheme.prefix)))


def test_oh2_charges():
    for scheme in iter_oh2_sto3g_gaussian_schemes():
        yield check_oh2_charges, scheme

def check_oh2_charges(scheme):
    expected = {
        "hirsh": numpy.array([-0.11247028, 0.11247028]),
        "hirshi": numpy.array([-0.17877806, 0.17877806]),
        "isa": numpy.array([-0.19463722, 0.19463722]),
    }
    scheme.do_charges()
    assert(abs(scheme.charges - expected[scheme.prefix]).max() < 1e-3)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_charges.txt" % scheme.prefix)))


def test_h_charges():
    for scheme in iter_h_sto3g_gaussian_schemes():
        yield check_h_charges, scheme

def check_h_charges(scheme):
    scheme.do_charges()
    assert(len(scheme.charges)==1)
    assert(abs(scheme.charges[0]) < 1e-3)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_charges.txt" % scheme.prefix)))


def test_hf_spin_charges():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_hf_spin_charges, scheme

def check_hf_spin_charges(scheme):
    scheme.do_spin_charges()
    assert(abs(scheme.spin_charges).max() < 1e-4)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_spin_charges.txt" % scheme.prefix)))


def test_oh1_spin_charges():
    for scheme in iter_oh1_sto3g_gaussian_schemes():
        yield check_oh1_spin_charges, scheme

def check_oh1_spin_charges(scheme):
    expected = {
        "hirsh": numpy.array([0.96630257, 0.03366921]),
        "hirshi": numpy.array([0.97163764, 0.0283441]),
        "isa": numpy.array([0.97359927, 0.02639019]),
    }
    scheme.do_spin_charges()
    assert(abs(scheme.spin_charges.sum() - 1.0) < 1e-3)
    assert(abs(scheme.spin_charges - expected[scheme.prefix]).max() < 1e-3)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_spin_charges.txt" % scheme.prefix)))


def test_oh2_spin_charges():
    for scheme in iter_oh2_sto3g_gaussian_schemes():
        yield check_oh2_spin_charges, scheme

def check_oh2_spin_charges(scheme):
    expected = {
        "hirsh": numpy.array([1.00267536, -0.00270729]),
        "hirshi": numpy.array([1.00719451, -0.00721775]),
        "isa": numpy.array([1.00897117, -0.00899064]),
    }
    scheme.do_spin_charges()
    assert(abs(scheme.spin_charges.sum() - 1.0) < 1e-3)
    assert(abs(scheme.spin_charges - expected[scheme.prefix]).max() < 1e-3)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_spin_charges.txt" % scheme.prefix)))


def test_hf_dipoles():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_hf_dipoles, scheme

def check_hf_dipoles(scheme):
    expected = {
        "hirsh": numpy.array([
            [ 7.33577068e-06,  4.21284970e-06, -1.30711487e-01],
            [-1.08521488e-05, -2.42288956e-05, -8.07245252e-02],
        ]),
        "hirshi": numpy.array([
            [-8.37620473e-06, -1.14026188e-05, -6.96008618e-02],
            [-6.05490107e-06, -1.02980542e-05, -3.37416403e-02],
        ]),
        "isa": numpy.array([
            [-4.34535236e-06,  1.78233520e-05, -5.70034272e-02],
            [-1.13215867e-05, -5.19636426e-06, -1.97794421e-02],
        ]),
        "becke": numpy.array([
            [-3.68237416e-05, -8.65399822e-06, -1.39946633e-01],
            [ 1.36355238e-06,  2.55892227e-06,  4.77311556e-02],
        ]),
    }
    scheme.do_dipoles()
    assert(abs(scheme.dipoles - expected[scheme.prefix]).max() < 1e-3)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_dipoles.txt" % scheme.prefix)))


def test_hf_noble_radii():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_hf_noble_radii, scheme

def check_hf_noble_radii(scheme):
    expected = numpy.array([0.32871748, 0.2])
    scheme.do_noble_radii()
    assert(abs(scheme.noble_radii - expected).max() < 1e-3)


def test_hf_bond_orders():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_hf_bond_orders, scheme

def check_hf_bond_orders(scheme):
    expected_bond_orders = {
        "hirsh": numpy.array([
            [0.0, 1.17375318],
            [1.17375318, 0.0],
        ]),
        "hirshi": numpy.array([
            [0.0, 1.11411172],
            [1.11411172, 0.0],
        ]),
        "isa": numpy.array([
            [0.0, 1.09701044],
            [1.09701044, 0.0],
        ]),
        "becke": numpy.array([
            [0.0, 1.04323287],
            [1.04323287, 0.0],
        ]),
    }
    expected_valences = {
        "hirsh": numpy.array([1.17377819, 1.17399224]),
        "hirshi": numpy.array([1.11407561, 1.11407691]),
        "isa": numpy.array([1.09685209, 1.09681649]),
        "becke": numpy.array([1.0432785, 1.04326027]),
    }
    scheme.do_bond_orders()
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_bond_orders.txt" % scheme.prefix)))
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_valences.txt" % scheme.prefix)))
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_free_valences.txt" % scheme.prefix)))
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_overlap_matrices.txt" % scheme.prefix)))
    assert(abs(scheme.bond_orders.sum(axis=0) - scheme.valences).max() < 1e-2)
    assert(abs(scheme.bond_orders - expected_bond_orders[scheme.prefix]).max() < 1e-3)
    assert(abs(scheme.valences - expected_valences[scheme.prefix]).max() < 1e-3)
    assert(abs(scheme.free_valences).max() < 1e-3)


def test_oh1_bond_orders():
    for scheme in iter_oh1_sto3g_gaussian_schemes():
        yield check_oh1_bond_orders, scheme

def check_oh1_bond_orders(scheme):
    expected_bond_orders = {
        "hirsh": numpy.array([
            [0.0, 1.20357067],
            [1.20357067, 0.0],
        ]),
        "hirshi": numpy.array([
            [0.0, 1.13842368],
            [1.13842368, 0.0],
        ]),
        "isa": numpy.array([
            [0.0, 1.11859939],
            [1.11859939, 0.0],
        ]),
    }
    expected_valences = {
        "hirsh": numpy.array([2.13744632, 1.20467815]),
        "hirshi": numpy.array([2.08258865, 1.13920176]),
        "isa": numpy.array([[2.06647643, 1.1193233]]),
    }
    expected_free_valences = {
        "hirsh": numpy.array([0.93387566, 0.00110749]),
        "hirshi": numpy.array([9.44164968e-01, 7.78083661e-04]),
        "isa": numpy.array([[9.47877044e-01, 7.23917974e-04]]),
    }
    scheme.do_bond_orders()
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_bond_orders.txt" % scheme.prefix)))
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_valences.txt" % scheme.prefix)))
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_free_valences.txt" % scheme.prefix)))
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_overlap_matrices.txt" % scheme.prefix)))
    assert(abs(scheme.bond_orders - expected_bond_orders[scheme.prefix]).max() < 1e-3)
    assert(abs(scheme.valences - expected_valences[scheme.prefix]).max() < 1e-3)
    assert(abs(scheme.free_valences - expected_free_valences[scheme.prefix]).max() < 1e-3)


def test_oh2_bond_orders():
    for scheme in iter_oh2_sto3g_gaussian_schemes():
        yield check_oh2_bond_orders, scheme

def check_oh2_bond_orders(scheme):
    expected_bond_orders = {
        "hirsh": numpy.array([
            [0.0, 1.20223341],
            [1.20223341, 0.0],
        ]),
        "hirshi": numpy.array([
            [0.0, 1.13753611],
            [1.13753611, 0.0],
        ]),
        "isa": numpy.array([
            [0.0, 1.1172814],
            [1.1172814, 0.0],
        ]),
    }
    expected_valences = {
        "hirsh": numpy.array([2.13860969, 1.20513964]),
        "hirshi": numpy.array([2.08410701, 1.13998475]),
        "isa": numpy.array([[2.06796024, 1.11955637]]),
    }
    expected_free_valences = {
        "hirsh": numpy.array([0.93637628, 0.00290623]),
        "hirshi": numpy.array([0.9465709, 0.00244864]),
        "isa": numpy.array([[0.95067885, 0.00227497]]),
    }
    scheme.do_bond_orders()
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_bond_orders.txt" % scheme.prefix)))
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_valences.txt" % scheme.prefix)))
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_free_valences.txt" % scheme.prefix)))
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_overlap_matrices.txt" % scheme.prefix)))
    assert(abs(scheme.bond_orders - expected_bond_orders[scheme.prefix]).max() < 1e-3)
    assert(abs(scheme.valences - expected_valences[scheme.prefix]).max() < 1e-3)
    assert(abs(scheme.free_valences - expected_free_valences[scheme.prefix]).max() < 1e-3)


def test_hf_net_overlap():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_hf_net_overlap, scheme

def check_hf_net_overlap(scheme):
    expected = {
        "hirsh": numpy.array([
            [8.89564521, 0.24206497],
            [0.24206497, 0.62019607],
        ]),
        "hirshi": numpy.array([
            [8.96906633, 0.22546657],
            [0.22546657, 0.58001127],
        ]),
        "isa": numpy.array([
            [8.98880362, 0.21955338],
            [0.21955338, 0.57209327],
        ]),
        "becke": numpy.array([
            [9.06191354, 0.13850754],
            [0.13850754, 0.66111818],
        ]),

    }
    scheme.do_net_overlap()
    assert(abs(scheme.net_overlap - expected[scheme.prefix]).max() < 1e-2)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_net_overlap.txt" % scheme.prefix)))


def test_hf_mol_esp_cost():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_hf_mol_esp_cost, scheme

def check_hf_mol_esp_cost(scheme):
    expected_A = numpy.array([
        [6.48538806e-03, 6.01517359e-03, 1.84842712e-06, -1.13002199e-06,
         2.42079860e-04, 5.24182736e-07, -1.97722198e-06, 2.45763010e-04],
        [6.01517359e-03, 5.77420167e-03, 1.27717341e-06, -1.42144727e-06,
         1.13216795e-04, 2.64573017e-07, -2.28072390e-06, 1.34829027e-04],
        [1.84842712e-06, 1.27717341e-06, 8.35583232e-05, -1.89864116e-08,
         3.45861341e-07, 6.34890875e-05, 1.20691190e-07, 2.37151713e-07],
        [-1.13002199e-06, -1.42144727e-06, -1.89864116e-08, 8.40184485e-05,
         2.72719526e-07, 1.20691190e-07, 6.39798830e-05, 4.68388644e-08],
        [2.42079860e-04, 1.13216795e-04, 3.45861341e-07, 2.72719526e-07,
         7.78488314e-05, 1.58151368e-07, 2.44930293e-07, 5.78852051e-05],
        [5.24182736e-07, 2.64573017e-07, 6.34890875e-05, 1.20691190e-07,
         1.58151368e-07, 5.95626741e-05, 1.65163832e-07, 1.11021715e-07],
        [-1.97722198e-06, -2.28072390e-06, 1.20691190e-07, 6.39798830e-05,
         2.44930293e-07, 1.65163832e-07, 5.99986010e-05, 7.79501731e-08],
        [2.45763010e-04, 1.34829027e-04, 2.37151713e-07, 4.68388644e-08,
         5.78852051e-05, 1.11021715e-07, 7.79501731e-08, 5.77542761e-05]
    ])
    expected_B = numpy.array([
        -1.19596516e-04, -6.29794625e-05, 2.83959175e-08, 1.04535738e-07,
        -3.11187068e-05, -4.36903490e-09, -3.49417636e-09, -2.78998269e-05
    ])
    expected_C = 1.4012993451e-05
    scheme.do_esp_costfunction()
    assert(abs(scheme.mol_esp_cost.A - expected_A).max() < 1e-4)
    assert(abs(scheme.mol_esp_cost.B - expected_B).max() < 1e-4)
    assert(abs(scheme.mol_esp_cost.C - expected_C) < 1e-6)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "mol_esp_cost.txt")))


def test_hf_esp_tst():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_hf_esp_tst, scheme

def check_hf_esp_tst(scheme):
    scheme.do_esp_test()
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_esp_test.txt" % scheme.prefix)))


def test_hf_multipoles():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_hf_multipoles, scheme

def check_hf_multipoles(scheme):
    expected = {
        "hirsh": numpy.array([
            [-1.37746685e-01, -1.30701945e-01, 8.57183561e-07, 7.68314867e-07,
             2.81712348e-01, 5.05192834e-06, 4.25054861e-06, 4.41685726e-06,
             -2.32785674e-06, 7.59746483e-02, -2.51590276e-05, -2.21909463e-05,
             8.64250844e-07, 3.16682962e-06, -1.83179635e-05, -2.38483096e-06,
             2.64473521e-03, 5.98338310e-05, 4.24018463e-05, -3.45317150e-05,
             -9.05323213e-06, 5.01964908e-05, 3.06125471e-06, 9.83346579e-06,
             -8.73856382e-06],
            [1.37715652e-01, -8.08417010e-02, -1.70767559e-05, -2.73748792e-05,
             2.06651352e-03, -4.09549331e-05, -4.87312106e-05, 2.87825995e-05,
             -3.36386615e-05, 4.38524213e-02, -2.24187353e-05, 2.54836745e-05,
             9.22727996e-05, -1.33213801e-04, 1.10245223e-04, -4.28431856e-05,
             5.81367149e-02, 1.90109554e-04, 2.77918275e-04, 2.74064764e-05,
             -2.80784452e-04, 1.19769537e-04, 2.05932568e-04, 7.28648232e-05,
             1.93034339e-04],
        ]),
        "hirshi": numpy.array([
            [-1.94559195e-01, -6.95851746e-02, -1.88934769e-06, -6.08044968e-08,
             2.30939669e-01, 1.21783212e-05, 4.74444275e-06, 5.91589084e-06,
             -2.51164912e-06, 1.04066396e-01, -4.02241198e-05, -1.82904631e-05,
             -6.86593312e-06, 4.06399329e-06, -1.89223567e-05, -5.48142174e-06,
             -3.31026457e-04, 8.62945335e-05, 2.64424019e-05, -4.73518444e-06,
             -1.35405964e-05, 5.28217193e-05, 1.87475087e-05, -4.46440733e-06,
             -3.53662742e-06],
            [1.94530691e-01, -3.37205189e-02, -1.67196366e-05, -2.62212007e-05,
             2.61952103e-02, -5.00926516e-05, -5.49219095e-05, 8.82139233e-06,
             -3.09346240e-05, 3.33508760e-02, -6.74593997e-05, -4.37163983e-06,
             7.61490105e-05, -8.01744565e-05, 9.73410454e-05, -7.55870053e-05,
             1.05742494e-02, 1.99078570e-04, 2.76203302e-04, 2.32135983e-04,
             -6.10307365e-05, 1.22845189e-04, 4.14160625e-05, 1.30258495e-04,
             1.04170285e-04],
        ]),
        "isa": numpy.array([
            [-2.08416652e-01, -5.70238910e-02, -9.10545057e-06, -8.07198380e-06,
              2.26655941e-01, 2.15556001e-05, 1.65951423e-05, 2.67054227e-06,
              1.54156087e-07, 9.56431115e-02, -4.89972907e-05, -2.26646719e-05,
              1.10866725e-05, -2.08233895e-05, 1.07486592e-05, 1.47592488e-05,
              2.10721172e-02, 1.27281926e-04, 4.33202049e-05, -7.69403127e-05,
              3.90957181e-05, -7.33942098e-05, -1.08665086e-05, -4.18051210e-05,
              1.92057601e-07],
            [2.08375963e-01, -1.98941611e-02, -2.22058937e-06, -3.33026585e-06,
             3.28958448e-02, -9.79062949e-07, -4.60282111e-05, 2.84774059e-05,
             -8.59080491e-05, 2.53203897e-02, -2.37650298e-05, -1.29318517e-04,
             5.71990474e-05, -2.77752875e-04, -8.48455245e-06, 7.71049887e-05,
             -1.79129679e-02, 1.09191539e-04, -1.93352643e-04, 1.09704115e-04,
             -1.92112594e-04, 5.10911851e-05, -6.29398231e-05, 4.30353724e-06,
             -8.27103047e-05],
        ]),
        "becke": numpy.array([
            [-2.00405499e-01, -1.39878557e-01, 7.97601430e-06, -1.66742368e-05,
             4.71749360e-01, -2.71808643e-05, 3.79557338e-05, -1.42330237e-05,
             -3.62462997e-05, -4.08015634e-01, 8.47184812e-05, -4.13451388e-05,
             2.40882895e-04, 3.30862111e-04, 3.88526966e-05, 1.99544087e-05,
             7.09491154e-01, -1.41657838e-04, 1.13581113e-04, -1.93087267e-03,
             -2.01981812e-03, -5.83155378e-05, -4.48755519e-05, 2.38971482e-04,
             1.48753657e-05],
            [2.00383290e-01, 4.77336964e-02, -1.07708896e-06, -2.01473416e-06,
             7.43597856e-02, -1.92734842e-06, -5.76092145e-06, 1.02021925e-06,
             2.60903874e-07, -2.53915134e-02, 4.01062845e-06, -1.24359990e-05,
             1.47981662e-07, 6.00618913e-06, 1.01071425e-06, -1.06564610e-06,
             -1.98158787e-02, 5.01886764e-05, -2.56465922e-05, -9.68329724e-06,
             3.10097861e-05, 8.52911267e-07, -4.86396262e-06, 6.55863937e-06,
             1.92435494e-06],
        ])
    }
    scheme.do_charges()
    scheme.do_dipoles()
    scheme.do_multipoles()
    assert(abs(scheme.charges - scheme.multipoles[:,0]).max() < 1e-3)
    assert(abs(scheme.dipoles[:,0] - scheme.multipoles[:,2]).max() < 1e-3)
    assert(abs(scheme.dipoles[:,1] - scheme.multipoles[:,3]).max() < 1e-3)
    assert(abs(scheme.dipoles[:,2] - scheme.multipoles[:,1]).max() < 1e-3)
    assert(abs(scheme.multipoles - expected[scheme.prefix]).max() < 1e-2)
    assert(os.path.isfile(os.path.join(scheme.context.outdir, "%s_multipoles.txt" % scheme.prefix)))


def test_hf_overlap_matrices():
    for scheme in iter_hf_sto3g_gaussian_schemes():
        yield check_hf_overlap_matrices, scheme

def check_hf_overlap_matrices(scheme):
    scheme.do_charges()
    scheme.do_atgrids_overlap_matrix()
    num_dof = scheme.context.wavefn.num_orbitals
    full = numpy.zeros((num_dof,num_dof), float)
    dmat_to_full(scheme.context.wavefn.density_matrix, full)
    molecule = scheme.context.wavefn.molecule
    for i in xrange(molecule.size):
        overlap = scheme.atgrids[i].overlap_matrix
        error = abs(overlap-overlap.transpose()).max()
        assert(error<1e-10)
        population = (full*overlap).sum()
        check_population = scheme.populations[i]
        error = abs(population-check_population)
        assert(error<1e-10)
