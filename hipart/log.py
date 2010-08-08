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
    def __init__(self, label, n, f=sys.stdout, indent=0, width=80):
        self.i = 0
        self.n = n
        self.f = f
        self.indent = indent
        self.width = width
        self.f.write("%s%s: " % (" "*indent, label))
        self.f.flush()

    def __call__(self, inc=1):
        for counter in xrange(inc):
            if self.i % self.width == 0:
                self.f.write("\n%s" % (" "*self.indent))
            if self.i % 10 == 0 or self.i == self.n:
                if self.i == 0:
                    self.f.write(" 0% .")
                else:
                    self.f.write(". %i%% " % ((100*self.i)/self.n))
            else:
                self.f.write(".")
            if self.i == self.n:
                self.f.write("\n")
            self.f.flush()
            self.i += 1


bright = '\033[1;33m'
reset = '\033[0m'

class Log(object):
    def __init__(self):
        self.stack = []

    level = property(lambda self: len(self.stack))

    def begin(self, s):
        print "%s%sBEGIN%s %s" % ("  "*self.level, bright, reset, s)
        self.stack.append(s)

    def __call__(self, s):
        print "%s%s" % ("  "*self.level, s)

    def end(self):
        s = self.stack.pop()
        print "%s%sEND%s %s" % ("  "*self.level, bright, reset, s)

    def pb(self, s, num):
        return ProgressBar(s, num, indent=self.level*2)


log = Log()
