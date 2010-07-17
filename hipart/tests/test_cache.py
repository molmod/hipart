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


from utils import FakeOptions, setup_hf_sto3g_gaussian

from hipart.cache import HirshfeldCache
from hipart.context import Context

import shutil


def test_compute_stockholder_weights():
    tmpdir, fn_fchk, fn_densities = setup_hf_sto3g_gaussian()
    options = FakeOptions(110, 50, 1, 1e-4, 500, True)
    context = Context(fn_fchk, options)
    cache = HirshfeldCache.new_from_args(context, [fn_densities])
    cache.do_atgrids()
    cache.do_proatomfns()
    h0 = cache._compute_atweights(cache.atgrids[0], 0)
    h1 = cache._compute_atweights(cache.atgrids[0], 1)
    assert(abs(h1+h0-1).max() < 1e-10)
    cache.do_charges()
    cache.do_dipoles()
    shutil.rmtree(tmpdir)
