from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from builtins import int
from builtins import range
from future import standard_library
standard_library.install_aliases()
from ase.data import covalent_radii as radii
import numpy as np
import spglib
from ase import Atoms


def get_neighbors(
        atoms,
        points=None,
        cutoff_matrix=None):
    """ Returns the neighboring atoms within a specified cutoff matrix
    criteria for requested points.

    Use of the cutoff matrix provides more fine-tuned control
    over the interaction parameters.

    Parameters:
      atoms: ASE atoms-object
        Atoms object to return.

      points: list (N,) or None
        Points to locate neighboring points of. If not provided, all
        atom points in the atoms-object will be used.

      cutoff_matrix: ndarray (96, 96) or None
        A matrix of interaction distances for the first 96 elements
        of the periodic table. These interactions are separated into
        individual i, j interactions. If None, defaults from ASE
        will be used.

    Returns:
      neighbors: dict
        Keys of each point and an array of each neighboring atoms index.
    """

    if cutoff_matrix is None:

        r = radii.copy()
        # TODO: develop an SVM to parameterize this for me
        # Will need reliable training data or an unsupervised approach
        metals = [
            [47, 1.1],  # Ag
            [79, 1.2],  # Au
            [29, 1.1],  # Cu
            [77, 1.1],  # Ir
            [46, 1.1],  # Pd
            [78, 1.2],  # Pt
            [45, 1.1]]  # Rh

        adsorbates = [
            [1, 1.0],
            [6, 1.0],
            [7, 1.0],
            [8, 1.0],
            [16, 1.0]]

        for i, f in metals:
            r[i] *= f
        for i, f in adsorbates:
            r[i] *= f
        cutoff_matrix = np.zeros((96, 96))
        for i in range(96):
            for j in range(96):
                cutoff_matrix[i][j] = r[i] + r[j]
                cutoff_matrix[j][i] = r[i] + r[j]

    rcmax = cutoff_matrix.max()

    an = atoms.get_atomic_numbers()
    positions = atoms.get_positions()
    pbc = atoms.get_pbc()
    cell = atoms.get_cell()

    icell = np.linalg.pinv(cell)
    scaled = np.dot(positions, icell)
    scaled0 = scaled.copy()

    N = []
    for i in range(3):
        if pbc[i]:
            scaled0[:, i] %= 1.0
            v = icell[:, i]
            h = 1 / np.sqrt(np.dot(v, v))
            n = int(2 * rcmax / h) + 1
        else:
            n = 0
        N.append(n)

    offsets = (scaled0 - scaled).round().astype(int)
    positions0 = atoms.positions + np.dot(offsets, cell)
    natoms = len(atoms)
    indices = np.arange(natoms)

    if points is None:
        points = indices

    cutoffs = np.zeros(natoms)
    neighbors = {a: np.empty(0, int) for a in points}
    for n1 in range(-N[0], N[0] + 1):
        for n2 in range(-N[1], N[1] + 1):
            for n3 in range(-N[2], N[2] + 1):

                displacement = np.dot((n1, n2, n3), cell)

                for a in points:

                    for b in range(natoms):
                        cutoffs[b] = cutoff_matrix[an[a]][an[b]]

                    d = positions0 + displacement - positions0[a]
                    i = indices[(d**2).sum(1) < (cutoffs)**2]

                    if a in i:
                        i = np.delete(i, np.where(i == a)[0])

                    neighbors[a] = np.concatenate(
                        (neighbors[a], i))

    return neighbors


def get_primitive_cell(
        atoms,
        tol=1e-8):
    """ ASE atoms-object interface with spglib primitive cell finder:
    https://atztogo.github.io/spglib/python-spglib.html#python-spglib

    Parameters:
      atoms: ASE atoms-object
        Atoms object to search for a primitive unit cell.

      tol: float
        Tolerance for floating point rounding errors.

    Returns:
      primitive cell: ASE atoms-object
        The primitive unit cell returned by spglib if one is found.
    """

    lattice = atoms.cell
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()

    cell = (lattice, positions, numbers)

    _lattice, _positions, _numbers = spglib.find_primitive(
        cell,
        symprec=tol)

    atoms = Atoms(symbols=_numbers, cell=_lattice)
    atoms.set_scaled_positions(_positions)

    return atoms


def get_symmetry(
        atoms,
        tol=1e-8):
    """ ASE atoms-object interface with spglib symmetry finder:
    https://atztogo.github.io/spglib/python-spglib.html#python-spglib

    Parameters:
      atoms: ASE atoms-object
        Atoms object to search for symmetric structures of.

      tol: float
        Tolerance for floating point rounding errors.

    Returns:
      symmetry operations: ndarray (N, N)
        Symmetry operations from spglib.
    """

    lattice = atoms.cell
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()

    cell = (lattice, positions, numbers)

    return spglib.get_symmetry(cell, symprec=tol)


def get_unique_coordinates(
        atoms,
        axis=2,
        direct=True,
        tol=1e-5):
    """ Return unique coordinate values of a given atoms object
    for a specified axis.

    Parameters:
      atoms: ASE atoms-object
        Atoms object to search for unique values along.

      axis: int
        Value of 0, 1, or 2 associated with x, y, and z coordinates.

      direct: bool
        Whether to use direct coordinates or Cartesian.

      tol: float
        The tolerance to search for unique values within.

    Returns:
      values: ndarray (N,)
        Array of unique values.
    """

    if direct:
        positions = atoms.get_scaled_positions()
    else:
        positions = atoms.positions

    values = [positions[0][axis]]
    for d in positions[1:, axis]:
        if not np.isclose(d, values, rtol=tol).any():
            values += [d]

    return np.array(values)


def vector_search(
        miller_index,
        limits):
    """ Try to find surface vectors within given range
        need to satisfy three conditions:
         - P.(h, k, l) = 0
         - Q.(h, k, l) = 0
         - PxQ = (h, k, l)

    Parameters:
      miller_index: list (3,)
        Three integers representing the miller indices.

      limits: list (2,)
        limits to search for surface vectors within.

    Returns:
      surface vectors: ndarray (2, 3)
        Array of surface vectors.
    """

    h, k, l = miller_index
    i0, i1 = limits

    for p1 in range(i0, i1 + 1):
        for p2 in range(i0, i1 + 1):
            for p3 in range(i0, i1 + 1):

                if(h * p1 + k * p2 + l * p3 != 0):
                    continue

                for q1 in range(i0, i1 + 1):
                    for q2 in range(i0, i1 + 1):
                        for q3 in range(i0, i1 + 1):

                            if h * q1 + k * q2 + l * q3 != 0:
                                continue
                            if p1 * q2 - p2 * q1 != l:
                                continue
                            if p3 * q1 - p1 * q3 != k:
                                continue
                            if p2 * q3 - p3 * q2 != h:
                                continue

                            return np.array([[p1, p2, p3], [q1, q2, q3]])
