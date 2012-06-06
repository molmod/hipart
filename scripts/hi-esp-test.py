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


from hipart.opts import parse_command_line

usage = """%prog test charges and dipoles on ESP grid

These atomic charges and dipoles are computed like in the hi-charges.py and
hi-dipoles.py program. Then a ESP costfunction is constructed and it is tested
how well the partitioned charges and dipoles reproduce the ESP."""


context, scheme = parse_command_line(usage)
scheme.do_esp_test()
context.work.clean()
