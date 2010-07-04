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


from hipart.atoms import AtomTable
from hipart.lebedev_laikov import get_grid
from hipart.wavefn import load_wavefunction

import os, shutil


__all__ = ["ContextError", "Context"]


class ContextError(Exception):
    pass


class Context(object):
    """This class is an extension to the concept of global variables.

       It ensures that an entire context of data structures, files and
       directories is present and sane. One can pass around a context variable
       instead of the individual global variables.
    """

    version = 2

    def __init__(self, filename, options):
        self.wavefn = load_wavefunction(filename)
        self.options = options

        self.outdir = "%s.hipart" % self.wavefn.prefix
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        self.workdir = os.path.join(self.outdir, "work")
        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

        self.lebedev_xyz, self.lebedev_weights = get_grid(options.lebedev)

    num_lebedev = property(lambda self: len(self.lebedev_weights))

    def check_tag(self, rs):
        """Make sure our context is compatible with the data in the workdir."""

        self.tag = "contextversion=%i lebedev=%i mol_lebedev=%i r_low=%.2e r_high=%.2e r_steps=%i filename=%s" % (
            self.version, self.num_lebedev, self.options.mol_lebedev,
            rs.min(), rs.max(), len(rs), os.path.basename(self.wavefn.filename),
        )

        context_fn = os.path.join(self.workdir, "context")
        if os.path.isfile(context_fn):
            f = file(context_fn)
            existing_tag = ("".join(f)).strip()
            f.close()
            if existing_tag != self.tag:
                raise ContextError("The existing work directory contains incompatible data. Trash it!")
        else:
            f = file(context_fn, "w")
            print >> f, self.tag
            f.close()

    def clean(self):
        if self.options.clean >= 3:
            print "Removing the entire work directory."
            shutil.rmtree(self.workdir)
            return
        if self.options.clean >= 2:
            print "Cleaning up binary files in work directory."
            os.system("rm -f %s" % os.path.join(self.workdir, "*.bin"))
        if self.options.clean >= 1:
            print "Cleaning up text files in work directory."
            os.system("rm -f %s" % os.path.join(self.workdir, "*.txt"))


