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
from hipart.log import log
from hipart.wavefn import load_wavefunction
from hipart.work import Work
from hipart.io import Output

import os


__all__ = ["ContextError", "Context", "Options"]


class Options(object):
    def __init__(self, lebedev=110, do_clean=False, do_work=True, do_output=True, verbose=True, threshold=1e-4, max_iter=500, fix_total_charge=True):
        self.lebedev = lebedev
        self.do_clean = do_clean
        self.do_work = do_work
        self.do_output = do_output
        self.verbose = verbose
        self.threshold = threshold
        self.max_iter = max_iter
        self.fix_total_charge = fix_total_charge


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
        log.verbose = options.verbose
        log.begin("Loading Electronic structure")
        self.wavefn = load_wavefunction(filename)
        self.wavefn.log()
        log.end()
        self.options = options

        outdir = "%s.hipart" % self.wavefn.prefix
        if options.do_output:
            self.output = Output(outdir, self.wavefn.molecule.numbers)
        else:
            self.output = Output()

        if options.do_work:
            workdir = os.path.join(outdir, "work")
        else:
            workdir = None
        self.work = Work(workdir, do_clean=options.do_clean)

    num_lebedev = property(lambda self: len(self.lebedev_weights))

    def check_tag(self, extra):
        """Make sure our context is compatible with the data in the workdir."""

        if not self.work.active:
            return

        tag_attributes = {
            "contextversion": "%i" % self.version,
            "filename": "%s" % os.path.basename(self.wavefn.filename),
        }
        tag_attributes.update(extra)

        context_fn = os.path.join(self.work.directory, "context")
        if os.path.isfile(context_fn):
            f = file(context_fn)
            existing_tag = ("".join(f)).strip()
            f.close()
            existing_attributes = dict(word.split("=") for word in existing_tag.split())
            for key, val in existing_attributes.iteritems():
                check = tag_attributes.get(key)
                if check is None:
                    tag_attributes[key] = val
                    continue
                if check != val:
                    message = [
                        "The existing work directory (%s) contains incompatible data.\n" % self.work.directory,
                        "Try using the --clean option once.\n",
                        "The following mismatch was detected in the work directory:\n",
                        "'%s' (found in work) versus '%s' (current script) for property '%s'" % (val, check, key),
                    ]
                    raise ContextError("".join(message))

        tag = " ".join("%s=%s" % (key, val) for key, val in sorted(tag_attributes.iteritems()))
        f = file(context_fn, "w")
        print >> f, tag
        f.close()
