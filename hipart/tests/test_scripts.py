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


import os, sys, subprocess

from hipart.tests.utils import setup_gaussian


# These tests just run the scripts to see if they do not crash on a simple
# example.


def test_scripts():
    schemes = ["becke", "hirsh", "hirshi", "isa"]
    fns_script = [
        "hi-bond-orders.py", "hi-charges.py", "hi-dipoles.py",
        "hi-esp-test.py", "hi-multipoles.py", "hi-net-overlap.py",
        "hi-overlap-matrices-orb.py", "hi-spin-charges.py",
    ]
    for scheme in schemes:
        for fn_script in fns_script:
            yield check_script, fn_script, scheme


def check_script(fn_script, scheme):
    fn_script = os.path.abspath(os.path.join("scripts", fn_script))
    tmpdir, fn_fchk, fn_densities = setup_gaussian("hf_sto3g", "sto3g")
    if scheme in ["hirsh", "hirshi"]:
        args = (fn_script, fn_fchk, scheme, fn_densities)
    else:
        args = (fn_script, fn_fchk, scheme)
    retcode = run(args)
    assert(retcode==0)


def run(args):
    f = file(args[0], "r")
    mod_args = ("/usr/bin/env", "python", "-") + args[1:]
    proc = subprocess.Popen(mod_args, stdin=f, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=os.getcwd())
    outdata, errdata = proc.communicate()
    f.close()
    print outdata
    print errdata
    return proc.returncode
