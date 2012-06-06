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


from hipart.schemes import scheme_classes
from hipart.context import Context, Options as FakeOptions

import os, tempfile, shutil, glob


__all__ = [
    "setup_gaussian", "iter_hf_sto3g_gaussian_schemes",
    "iter_hf_sto3g_gaussian_schemes_opts", "iter_oh1_sto3g_gaussian_schemes",
    "iter_oh2_sto3g_gaussian_schemes", "iter_h_sto3g_gaussian_schemes",
]


def clean_txt(directory):
    if directory is not None:
        for fn_txt in glob.glob(os.path.join(directory, "*.txt")):
            os.remove(fn_txt)


def setup_gaussian(fchk_name, densities_name):
    tmpdir = tempfile.mkdtemp("hipart")
    if not os.path.isdir("input"):
        raise IOError("Input directory with test files is not present")
    fn_fchk = os.path.join(tmpdir, "gaussian.fchk")
    shutil.copy("input/%s.fchk" % fchk_name, fn_fchk)
    fn_densities = os.path.join(tmpdir, "densities.txt")
    shutil.copy("input/densities_%s.txt" % densities_name, fn_densities)
    return tmpdir, fn_fchk, fn_densities

def iter_hf_sto3g_gaussian_schemes():
    options = FakeOptions()

    tmpdir, fn_fchk, fn_densities = setup_gaussian("hf_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    shutil.rmtree(tmpdir)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("hf_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['hirshi'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirshi'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    shutil.rmtree(tmpdir)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("hf_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['isa'].new_from_args(context, ["2e-5", "20.0", "100"])
    clean_txt(context.output.directory)
    yield scheme_classes['isa'].new_from_args(context, [])
    shutil.rmtree(tmpdir)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("hf_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['becke'].new_from_args(context, ["3", "2e-5", "20.0", "100"])
    clean_txt(context.output.directory)
    yield scheme_classes['becke'].new_from_args(context, ["3"])
    clean_txt(context.output.directory)
    yield scheme_classes['hirshi'].new_from_args(context, [fn_densities])
    shutil.rmtree(tmpdir)


def iter_hf_sto3g_gaussian_schemes_opts():
    options = FakeOptions(do_clean=True)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("hf_sto3g", "sto3g")
    for do_random in True, False:
        options.do_random = do_random
        for do_work in True, False:
            options.do_work = do_work
            for do_output in True, False:
                options.do_output = do_output
                for save_mem in True, False:
                    options.save_mem = save_mem
                    context = Context(fn_fchk, options)
                    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
                    clean_txt(context.output.directory)
    shutil.rmtree(tmpdir)


def iter_oh1_sto3g_gaussian_schemes():
    options = FakeOptions()

    tmpdir, fn_fchk, fn_densities = setup_gaussian("oh_rad1_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    shutil.rmtree(tmpdir)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("oh_rad1_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['hirshi'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirshi'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    shutil.rmtree(tmpdir)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("oh_rad1_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['isa'].new_from_args(context, ["2e-5", "20.0", "100"])
    clean_txt(context.output.directory)
    yield scheme_classes['isa'].new_from_args(context, [])
    shutil.rmtree(tmpdir)


def iter_oh2_sto3g_gaussian_schemes():
    options = FakeOptions()

    tmpdir, fn_fchk, fn_densities = setup_gaussian("oh_rad2_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    shutil.rmtree(tmpdir)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("oh_rad2_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['hirshi'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirshi'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    shutil.rmtree(tmpdir)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("oh_rad2_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['isa'].new_from_args(context, ["2e-5", "20.0", "100"])
    clean_txt(context.output.directory)
    yield scheme_classes['isa'].new_from_args(context, [])
    shutil.rmtree(tmpdir)


def iter_h_sto3g_gaussian_schemes():
    options = FakeOptions()

    tmpdir, fn_fchk, fn_densities = setup_gaussian("h_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    shutil.rmtree(tmpdir)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("h_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['hirshi'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirshi'].new_from_args(context, [fn_densities])
    clean_txt(context.output.directory)
    yield scheme_classes['hirsh'].new_from_args(context, [fn_densities])
    shutil.rmtree(tmpdir)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("h_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['isa'].new_from_args(context, ["2e-5", "20.0", "100"])
    clean_txt(context.output.directory)
    yield scheme_classes['isa'].new_from_args(context, [])
    shutil.rmtree(tmpdir)

    tmpdir, fn_fchk, fn_densities = setup_gaussian("h_sto3g", "sto3g")
    context = Context(fn_fchk, options)
    yield scheme_classes['becke'].new_from_args(context, ["2e-5", "20.0", "100"])
    clean_txt(context.output.directory)
    yield scheme_classes['becke'].new_from_args(context, [])
    shutil.rmtree(tmpdir)
