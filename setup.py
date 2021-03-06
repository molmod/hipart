#!/usr/bin/env python
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


import glob
from numpy.distutils.core import setup, Extension


setup(
    name='HiPart',
    version='0.004',
    description='HiPart is a tool to analyse molecular densities with the fuzzy atom partitioning',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    url='http://molmod.ugent.be/code/',
    packages = ['hipart', 'hipart.gint'],
    package_dir = {'hipart': 'hipart'},
    scripts=glob.glob("scripts/hi-*.py"),
    license = "GPLv3",
    ext_modules=[
        Extension(
            "hipart.llext",
            ["hipart/Lebedev-Laikov.F"]
        ),
        Extension(
            "hipart.ext",
            ["hipart/cubic_spline.c", "hipart/utils.c", "hipart/ext.pyf"]
        ),
        Extension(
            "hipart.gint.gint_ext",
            glob.glob("hipart/gint/*.c") + ["hipart/gint/gint_ext.pyf"],
            depends=glob.glob("hipart/gint/*.pyf.inc") +
                    glob.glob("hipart/gint/*.h") +
                    ["hipart/gint/gaux.inc"]
        ),
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Topic :: Science/Engineering :: Molecular Science'
    ],
)
