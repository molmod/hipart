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






import sys


__all__ = ["ProgressBar", "Log", "log"]


class ProgressBar(object):
    def __init__(self, label, n, indent=0, width=80, verbose=True):
        self.i = 0
        self.n = n
        self.indent = indent
        self.width = width
        self.verbose = verbose
        if self.verbose:
            sys.stdout.write("%s%s: " % (" "*indent, label))
            sys.stdout.flush()

    def __call__(self, inc=1):
        if not self.verbose:
            return
        for counter in xrange(inc):
            if self.i % self.width == 0:
                sys.stdout.write("\n%s" % (" "*self.indent))
            if self.i % 10 == 0 or self.i == self.n:
                if self.i == 0:
                    sys.stdout.write(" 0% .")
                else:
                    sys.stdout.write(". %i%% " % ((100*self.i)/self.n))
            else:
                sys.stdout.write(".")
            if self.i == self.n:
                sys.stdout.write("\n")
            sys.stdout.flush()
            self.i += 1


bright = '\033[1;33m'
reset = '\033[0m'

class Log(object):
    def __init__(self):
        self.stack = []
        self.verbose = True

    level = property(lambda self: len(self.stack))

    def begin(self, s):
        if self.verbose:
            print "%s%sBEGIN%s %s" % ("  "*self.level, bright, reset, s)
        self.stack.append(s)

    def __call__(self, s):
        if self.verbose:
            print "%s%s" % ("  "*self.level, s)

    def end(self):
        s = self.stack.pop()
        if self.verbose:
            print "%s%sEND%s %s" % ("  "*self.level, bright, reset, s)

    def pb(self, s, num):
        return ProgressBar(s, num, self.level*2, 80-self.level*2, self.verbose)


log = Log()
