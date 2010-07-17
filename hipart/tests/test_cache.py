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


from utils import iter_hf_sto3g_gaussian_caches


def test_compute_stockholder_weights():
    for cache in iter_hf_sto3g_gaussian_caches():
        yield check_compute_stockholder_weights, cache

def check_compute_stockholder_weights(cache):
    cache.do_atgrids()
    cache.do_proatomfns()
    h0 = cache._compute_atweights(cache.atgrids[0], 0)
    h1 = cache._compute_atweights(cache.atgrids[0], 1)
    error = abs(h1+h0-1).max()
    assert(error < 1e-10)
    cache.do_charges()
    cache.do_dipoles()
