# HiPart is a tool to analyse molecular densities with the hirshfeld partitioning scheme
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


from hipart.atoms import AtomTable
from hipart.cache import Cache
from hipart.tools import guess_density_type
from hipart.lebedev_laikov import get_grid

from molmod.io.gaussian03.fchk import FCHKFile

import os, shutil


__all__ = ["ContextError", "Context"]


class ContextError(Exception):
    pass


class Context(object):
    """This class is an extension to the concept of global variables.

    It ensures that an entire context of data structures, files and directories
    is present and sane. One can pass around a context variable instead of the
    individual variables.
    """

    version = 1

    def __init__(self, atom_table, fchk_fn, options):
        self.atom_table = atom_table
        self.fchk = FCHKFile(fchk_fn, field_labels=[
            "Charge", "Number of basis functions", "Dipole Moment"
        ])
        self.options = options

        prefix = fchk_fn.replace(".fchk", "")
        self.outdir = "%s.hipart" % prefix
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        self.workdir = os.path.join(self.outdir, "work")
        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

        if options.density is None:
            options.density = guess_density_type(self.fchk.lot)
        else:
            options.density = options.density

        self.lebedev_xyz, self.lebedev_weights = get_grid(options.lebedev)

        self.check_context()
        if (options.reference is None or os.path.samefile(options.reference, fchk_fn)):
            self.reference = None
        else:
            self.reference = Context(atom_table, options.reference, options)
            if self.tag != self.reference.tag:
                raise ContectError("The reference state must have the same context tag.")
        self.cache = Cache(self)

    num_lebedev = property(lambda self: len(self.lebedev_weights))

    def check_context(self):
        """Make sure our context is compatible with the data in the workdir."""

        self.tag = "contextversion=%i density=%s lebedev=%i mol_lebedev=%i" % (
            self.version, self.options.density, self.num_lebedev, self.options.mol_lebedev
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
            os.system("rm -f %s" % os.path.join(self.workdir, "*.cube"))


