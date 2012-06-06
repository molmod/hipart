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


import sys, time


__all__ = ["ProgressBar", "Log", "log"]


class ProgressBar(object):
    def __init__(self, label, n, indent=0, width=80, verbose=True):
        self.i = 0
        self.n = n
        self.indent = indent
        self.width = width
        self.verbose = verbose
        self.pos = 0
        if self.verbose:
            self.dump(label + ": ")

    def dump(self, s):
        if len(s) > self.width - self.pos:
            part1 = s[:self.width - self.pos]
            part2 = s[self.width - self.pos:]
            self.dump(part1)
            self.dump(part2)
            return
        if self.pos == 0:
            sys.stdout.write(" "*self.indent)
        sys.stdout.write(s)
        self.pos = (self.pos + len(s)) % self.width
        if self.pos == 0:
            sys.stdout.write("\n")
        sys.stdout.flush()

    def __call__(self, inc=1):
        if not self.verbose:
            return
        line = []
        for counter in xrange(inc):
            if self.i % 10 == 0 or self.i == self.n:
                if self.n == 0:
                    pct = 100
                else:
                    pct = (100*self.i)/self.n
                line.append("%i%%" % pct)
            else:
                line.append(".")
            if self.i == self.n:
                line.append("\n")
            self.i += 1
        self.dump("".join(line))


if sys.stdout.isatty():
    bright = '\033[1;33m'
    reset = '\033[0m'
else:
    bright = ''
    reset = ''


class Log(object):
    def __init__(self):
        self.stack = []
        self._verbose = False
        self.start = time.time()

    def set_verbose(self, verbose):
        if verbose and not self._verbose:
            print "---TIME--- "+"-"*33+"LOG"+"-"*33
        self._verbose = verbose

    level = property(lambda self: len(self.stack))

    def _section(self, key, s):
        space = 69-self.level*2
        line = "%s%s%s %s" % (bright, key, reset, s)
        prefix = "%10.2f" % (time.time() - self.start)
        start = 0
        rows = (len(line) - len(bright) - len(reset))/space+1
        for i in xrange(rows):
            end = start + space
            if i==0:
                end += len(bright) + len(reset)
            print "%s %s%s" % (prefix, "  "*self.level, line[start:end])
            prefix = "          "
            start = end

    def begin(self, s):
        if self._verbose:
            self._section("BEGIN", s)
        self.stack.append(s)

    def __call__(self, s):
        if self._verbose:
            now = time.time() - self.start
            print "%10.2f %s%s" % (now, "  "*self.level, s)

    def end(self):
        s = self.stack.pop()
        if self._verbose:
            self._section("END", s)

    def pb(self, s, num):
        indent = self.level*2+11
        return ProgressBar(s, num, indent, 80-indent, self._verbose)


log = Log()
