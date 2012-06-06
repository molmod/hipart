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


from hipart.atomdb import charge_to_label, label_to_charge, parse_numbers

def test_charge_to_label():
    assert(charge_to_label(0) == "neut")
    assert(charge_to_label(2) == "pos2")
    assert(charge_to_label(-3) == "neg3")

def test_label_to_charge():
    assert(label_to_charge("neut") == 0)
    assert(label_to_charge("pos1") == 1)
    assert(label_to_charge("neg2") == -2)

def test_parse_numbers():
    assert(parse_numbers("1-3")==[1,2,3])
    assert(parse_numbers("5,7,8")==[5,7,8])
    assert(parse_numbers("1,6-9")==[1,6,7,8,9])
