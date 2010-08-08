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


from hipart.log import log

from molmod.periodic import periodic

import numpy, os


__all__ = [
    "Output",
    "dump_atom_scalars", "load_atom_scalars",
    "dump_atom_vectors", "load_atom_vectors",
    "dump_atom_matrix", "load_atom_matrix",
    "dump_atom_fields", "load_atom_fields",
    "dump_overlap_matrices",
]


class Output(object):
    def __init__(self, directory=None, numbers=None):
        self.directory = directory
        self.numbers = numbers
        if self.active and not os.path.isdir(directory):
            os.makedirs(directory)

    active = property(lambda self: self.directory is not None)

    def dump_atom_scalars(self, filename, scalars, name):
        if self.active:
            filename = os.path.join(self.directory, filename)
            dump_atom_scalars(filename, scalars, name, self.numbers)
            log("Written %s" % filename)

    def dump_atom_vectors(self, filename, vectors, name):
        if self.active:
            filename = os.path.join(self.directory, filename)
            dump_atom_vectors(filename, vectors, name, self.numbers)
            log("Written %s" % filename)

    def dump_atom_matrix(self, filename, matrix, name):
        if self.active:
            filename = os.path.join(self.directory, filename)
            dump_atom_matrix(filename, matrix, name, self.numbers)
            log("Written %s" % filename)

    def dump_atom_fields(self, filename, table, labels, name):
        if self.active:
            filename = os.path.join(self.directory, filename)
            dump_atom_fields(filename, table, labels, name, self.numbers)
            log("Written %s" % filename)

    def dump_overlap_matrices(self, filename, overlap_matrices):
        if self.active:
            filename = os.path.join(self.directory, filename)
            dump_overlap_matrices(filename, overlap_matrices, self.numbers)
            log("Written %s" % filename)

    def dump_esp_test(self, filename, dipole_q, dipole_p, dipole_qp, dipole_qm, mol_esp_cost, charges, dipoles):
        if self.active:
            filename = os.path.join(self.directory, filename)
            dump_esp_test(filename, dipole_q, dipole_p, dipole_qp, dipole_qm, mol_esp_cost, charges, dipoles)
            log("Written %s" % filename)

    def dump_esp_cost(self, filename, esp_cost):
        if self.active:
            filename = os.path.join(self.directory, filename)
            esp_cost.write_to_file(filename)
            log("Written %s" % filename)


def _iter_symbols_numbers(numbers, N):
    """Iterate over the given atom numbers and give additional info.

       Arguments:
        | ``numbers``  --  An array with atom numbers or None.
        | ``N``  --  The total number of atoms. (This is not derived from the
                     array with atom numbers as it may be None.)
    """
    for i in xrange(N):
        if numbers is None:
            number = 0
            symbol = "?"
        else:
            number = numbers[i]
            symbol = periodic[number].symbol
        yield i, symbol, number


def dump_atom_scalars(filename, scalars, name, numbers=None):
    """Dump an array of scalar atomic quantities into a text file.

       Arguments:
        | ``filename``  --  The file to dump in.
        | ``scalars``  --  An array of scalars, e.g. atomic charges.

       Optional arguments:
        | ``numbers``  --  An array with atomic numbers to decorate the file
        | ``name``  --  the name of the quantity to decorate the file
    """
    f = file(filename, "w")
    print >> f, "number of atoms: %i" % len(scalars)
    print >> f, "  i        Z  %s" % name.rjust(10)
    print >> f, "--------------------------------"
    for i, symbol, number in _iter_symbols_numbers(numbers, len(scalars)):
        print >> f, "% 3i  %2s  % 3i   % 15.12f" % (
            i, symbol, number, scalars[i]
        )
    print >> f, "--------------------------------"


def load_atom_scalars(filename):
    """Load atomic scalars written with :func:`dump_atom_scalars`.

       Argument:
        | ``filename``  --  The file to load from.

       Returns the array with scalars.
    """
    f = file(filename)
    # read the number of atoms
    line = f.next()
    N = int(line[line.rfind(" "):])
    # read the scalars
    scalars = numpy.zeros(N, float)
    f.next()
    f.next()
    for i in xrange(N):
        line = f.next()
        words = line.split()
        scalars[i] = float(words[3])
    f.close()
    return scalars


def dump_atom_vectors(filename, vectors, name, numbers=None):
    """Dump an array of atomic 3D-vector quantities into a text file.

       Arguments:
        | ``filename``  --  The file to dump in.
        | ``vectors``  --  An array of 3D-vectors, e.g. atomic dipoles.

       Optional arguments:
        | ``numbers``  --  An array with atomic numbers to decorate the file
        | ``name``  --  the name of the quantity to decorate the file
    """
    names = (
        (name+"-X").center(15),
        (name+"-Y").center(15),
        (name+"-Z").center(15),
        (name+"-norm").center(15),
    )
    f = file(filename, "w")
    print >> f, "number of atoms:", len(vectors)
    print >> f, "  i        Z  %s %s %s %s" % names
    print >> f, "-------------------------------------------------------------------------------"
    for i, symbol, number in _iter_symbols_numbers(numbers, len(vectors)):
        print >> f, "% 3i  %2s  % 3i  % 15.12f % 15.12f % 15.12f % 15.12f" % (
            i, symbol, number, vectors[i,0], vectors[i,1],
            vectors[i,2], numpy.linalg.norm(vectors[i]),
        )
    f.close()


def load_atom_vectors(filename):
    """Load atomic vectors written with :func:`dump_atom_vectors`.

       Argument:
        | ``filename``  --  The file to load from.

       Returns the array with vectors.
    """
    f = file(filename)
    # read the number of atoms
    line = f.next()
    N = int(line[line.rfind(" "):])
    # read the vectors
    vectors = numpy.zeros((N,3), float)
    f.next()
    f.next()
    for i in xrange(N):
        line = f.next()
        words = line.split()
        vectors[i, 0] = float(words[3])
        vectors[i, 1] = float(words[4])
        vectors[i, 2] = float(words[5])
    f.close()
    return vectors


def dump_atom_matrix(filename, matrix, name, numbers=None):
    """Dump a 2D-array of atomic pair quantities into a text file.

       Arguments:
        | ``filename``  --  The file to dump in.
        | ``matrix``  --  A 2D-array of atomic pair quantities, e.g. bond
                          orders.

       Optional arguments:
        | ``numbers``  --  An array with atomic numbers to decorate the file
        | ``name``  --  the name of the quantity to decorate the file
    """
    name = name.center(15)
    f = file(filename, "w")
    print >> f, "number of atoms:", len(matrix)
    print >> f, "%s | %s" % (name, " ".join(
        "    % 3i %2s     " % (key[0], key[1]) for key
        in _iter_symbols_numbers(numbers, len(matrix))
    ))
    print >> f, "----------------+-"+"-"*(1+16*len(matrix))
    for i, symbol, number in _iter_symbols_numbers(numbers, len(matrix)):
        print >> f, "% 3i  %2s  % 3i    | %s" % (
            i, symbol, number, " ".join("% 15.12f" % val for val in matrix[i])
        )
    f.close()


def load_atom_matrix(filename):
    """Load atomic pair quantities written with :func:`dump_atom_matrix`.

       Argument:
        | ``filename``  --  The file to load from.

       Returns the square array with pair quantities.
    """
    f = file(filename)
    # read the number of atoms
    line = f.next()
    N = int(line[line.rfind(" "):])
    # read the matrix
    matrix = numpy.zeros((N,N), float)
    f.next()
    f.next()
    for i in xrange(N):
        line = f.next()
        words = line.split()[4:]
        for j in xrange(N):
            matrix[i,j] = float(words[j])
    f.close()
    return matrix


def dump_atom_fields(filename, table, labels, name, numbers=None):
    """Dump a table with multiple scalar atomic quantities into a text file.

       Arguments:
        | ``filename``  --  The file to dump in.
        | ``table``  --  A 2D-array of atomic scalar quantities, e.g. the atomic
                         multipole expansions. Each column corresponds to a
                         quantity and each row corresponds to an atom.
        | ``labels``  --  The labels of the atomic quantities. The number of
                          labels and the number of columns in the table must be
                          the same.

       Optional arguments:
        | ``numbers``  --  An array with atomic numbers to decorate the file
        | ``name``  --  the name of the quantity to decorate the file
    """
    if len(labels) != table.shape[1]:
        raise ValueError("The number of labels must be equal to the number of atoms in the table.")
    name = name.center(15)
    f = file(filename, "w")
    print >> f, "number of atoms:", len(table)
    print >> f, "number of fields:", len(labels)
    print >> f, "%s | %s" % (name, " ".join(label.center(15) for label in labels))
    print >> f, "----------------+-"+"-"*(1+16*len(labels))
    for i, symbol, number in _iter_symbols_numbers(numbers, len(table)):
        print >> f, "% 3i  %2s  % 3i    | %s" % (
            i, symbol, number, " ".join("% 15.12f" % val for val in table[i])
        )
    f.close()


def load_atom_fields(filename):
    """Load multiple atomic quantities written with :func:`dump_atom_fields`.

       Argument:
        | ``filename``  --  The file to load from.

       Returns a tuple with two results: (i) the table with atomic scalar
       quantities and (ii) the labels from the table header.
    """
    f = file(filename)
    # read the number of atoms
    line = f.next()
    N = int(line[line.rfind(" "):])
    # read the labels
    f.next()
    line = f.next()
    line = line[line.find("|")+1:]
    labels = line.split()
    # read the table
    table = numpy.zeros((N,len(labels)), float)
    f.next()
    for i in xrange(N):
        line = f.next()
        words = line.split()[4:]
        for j, word in enumerate(words):
            table[i,j] = float(word)
    f.close()
    return table, labels


def dump_overlap_matrices(filename, overlap_matrices, numbers=None):
    """Dump a the atomic overlap matrices into a text file.

       Arguments:
        | ``filename``  --  The file to dump in.
        | ``overlap_matrices``  --  A list with (square) atomic overlap
                                    matrices.

       Optional arguments:
        | ``numbers``  --  An array with atomic numbers to decorate the file
    """
    f = file(filename, "w")
    print >> f, "number of orbitals:", len(overlap_matrices[0])
    print >> f, "number of atoms: ", len(overlap_matrices)
    for i, symbol, number in _iter_symbols_numbers(numbers, len(overlap_matrices)):
        print >> f, "Atom % 3i  %2s  % 3i" % (i, symbol, number)
        matrix = overlap_matrices[i]
        for row in matrix:
            print >> f, " ".join("% 15.10e" % value for value in row)
    f.close()


def dump_esp_test(filename, dipole_q, dipole_p, dipole_qp, dipole_qm, mol_esp_cost, charges, dipoles):
    f = file(filename, "w")
    print >> f, "Reproduction of the molecular dipole"
    print >> f, "-------------------------------------------------------------------------------"
    print >> f, "                  Dipole-X        Dipole-Y        Dipole-Z       Dipole-norm"
    print >> f, "-------------------------------------------------------------------------------"
    print >> f, "charges (q)   % 15.12f % 15.12f % 15.12f % 15.12f" % (
        dipole_q[0], dipole_q[1], dipole_q[2], numpy.linalg.norm(dipole_q),
    )
    print >> f, "dipoles (p)   % 15.12f % 15.12f % 15.12f % 15.12f" % (
        dipole_p[0], dipole_p[1], dipole_p[2], numpy.linalg.norm(dipole_p),
    )
    print >> f, "q and p       % 15.12f % 15.12f % 15.12f % 15.12f" % (
        dipole_qp[0], dipole_qp[1], dipole_qp[2], numpy.linalg.norm(dipole_qp),
    )
    print >> f, "total density % 15.12f % 15.12f % 15.12f % 15.12f" % (
        dipole_qm[0], dipole_qm[1], dipole_qm[2], numpy.linalg.norm(dipole_qm),
    )
    print >> f, "-------------------------------------------------------------------------------"
    print >> f
    print >> f, "Reproduction of the external molecular ESP"
    print >> f, "-------------------------------------------------------------"
    print >> f, "                     RMSD             RMS       CORRELATION"
    print >> f, "-------------------------------------------------------------"
    print >> f, "charges (q)      % 10.5e    % 10.5e      % 5.2f" % (
        mol_esp_cost.rmsd(charges),
        mol_esp_cost.model_rms(charges),
        mol_esp_cost.correlation(charges),
    )
    print >> f, "dipoles (p)      % 10.5e    % 10.5e      % 5.2f" % (
        mol_esp_cost.rmsd(None, dipoles),
        mol_esp_cost.model_rms(None, dipoles),
        mol_esp_cost.correlation(None, dipoles),
    )
    print >> f, "q and p          % 10.5e    % 10.5e      % 5.2f" % (
        mol_esp_cost.rmsd(charges, dipoles),
        mol_esp_cost.model_rms(charges, dipoles),
        mol_esp_cost.correlation(charges, dipoles),
    )
    print >> f, "total density                    % 10.5e" % mol_esp_cost.rms
    print >> f, "-------------------------------------------------------------"
    f.close()
    log("Written %s" % filename)
