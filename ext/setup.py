# HiPart is a software toolkit to analyse molecular densities with the hirshfeld partitioning scheme.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension

setup(
    name='HiPart',
    version='0.001',
    description='HiPart is a tool to analyse molecular densities with the hirshfeld partitioning scheme',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    url='http://molmod.ugent.be/code/',
    ext_modules=[
        Extension("hipart.lebedev_laikov_ext", ["lib/Lebedev-Laikov.F"]),
        Extension("hipart.cubic_spline_ext", ["lib/cubic_spline_ext.c","lib/cubic_spline_ext.pyf"]),
    ],
    packages=['hipart'],
    package_dir={'hipart': 'lib'},
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

