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

    def __call__(self):
        if self.i % self.width == 0:
            self.f.write("\n%s" % (" "*self.indent))
        if self.i % 10 == 0 or self.i == self.n:
            self.f.write(" %i%% " % ((100*self.i)/self.n))
        else:
            self.f.write(".")
        if self.i == self.n:
            self.f.write("\n")
        self.f.flush()
        self.i += 1


class Log(object):
    def __init__(self):
        self.level = 0

    def begin(self, s):
        print "%sBEGIN %s" % (" "*self.level, s)
        self.level += 1

    def __call__(self, s):
        print "%s%s" % (" "*self.level, s)

    def end(self, s):
        self.level -= 1
        if self.level < 0: self.level = 0
        print "%sEND %s" % (" "*self.level, s)

    def pb(self, s, num):
        return ProgressBar(s, num, indent=self.level)


log = Log()


