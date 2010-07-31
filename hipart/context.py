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
        self.wavefn.log()
        self.options = options

        self.outdir = "%s.hipart" % self.wavefn.prefix
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        self.workdir = os.path.join(self.outdir, "work")
        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

        self.lebedev_xyz, self.lebedev_weights = get_grid(options.lebedev)

    num_lebedev = property(lambda self: len(self.lebedev_weights))

    def check_tag(self, extra):
        """Make sure our context is compatible with the data in the workdir."""

        tag_attributes = {
            "contextversion": "%i" % self.version,
            "lebedev": "%i" % self.num_lebedev,
            "filename": "%s" % os.path.basename(self.wavefn.filename),
        }
        tag_attributes.update(extra)
        if hasattr(self.options, "mol_lebedev"):
            tag_attributes["mol_lebedev"] = "%i" % self.options.mol_lebedev

        context_fn = os.path.join(self.workdir, "context")
        if os.path.isfile(context_fn):
            f = file(context_fn)
            existing_tag = ("".join(f)).strip()
            f.close()
            existing_attributes = dict(word.split("=") for word in existing_tag.split())
            for key, val in existing_attributes.iteritems():
                check = tag_attributes.get(key)
                if check is None:
                    continue
                if check != val:
                    raise ContextError("The existing work directory contains incompatible data. Trash it!")

        tag = " ".join("%s=%s" % (key, val) for key, val in sorted(tag_attributes.iteritems()))
        f = file(context_fn, "w")
        print >> f, tag
        f.close()

    def clean(self):
        if self.options.clean:
            shutil.rmtree(self.workdir)
