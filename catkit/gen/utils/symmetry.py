from catkit import Gratoms
import numpy as np
import spglib


def get_spglib_cell(atoms, primitive=False, idealize=True, tol=1e-5):
    """Atoms object interface with spglib primitive cell finder:
    https://atztogo.github.io/spglib/python-spglib.html#python-spglib

    Parameters
    ----------
    atoms : object
        Atoms object to search for a primitive unit cell.
    primitive : bool
        Reduce the atoms object into a primitive form.
    idealize : bool
        Convert the cell into the spglib standardized form.
    tol : float
        Tolerance for floating point rounding errors.

    Returns
    -------
    primitive cell : object
        The primitive unit cell returned by spglib if one is found.
    """
    lattice = atoms.cell
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()

    cell = (lattice, positions, numbers)
    cell = spglib.standardize_cell(
        cell, to_primitive=primitive, no_idealize=~idealize, symprec=tol)

    if cell is None:
        return atoms

    _lattice, _positions, _numbers = cell
    atoms = Gratoms(symbols=_numbers, cell=_lattice, pbc=atoms.pbc)
    atoms.set_scaled_positions(_positions)

    return atoms


def get_point_group(atoms, tol=1e-8):
    """Return the point group operations of a systems."""
    rotations, translations = get_symmetry(atoms, tol=tol)
    point_group = spglib.get_pointgroup(rotations)[0].strip()

    laue = ['-1', '2/m', 'mmm', '4/m', '4/mmm',
            '-3', '-3m', '6/m', '6/mmm', 'm-3', 'm-3m']

    is_laue = point_group in laue

    return point_group, is_laue


def get_symmetry(atoms, tol=1e-8):
    """Atoms object interface with spglib symmetry finder:
    https://atztogo.github.io/spglib/python-spglib.html#python-spglib

    Parameters
    ----------
    atoms : object
        Atoms object to search for symmetric structures of.
    tol : float
        Tolerance for floating point rounding errors.

    Returns
    -------
    symmetry operations: ndarray (n, n)
        Symmetry operations from spglib.
    """
    lattice = atoms.cell
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()

    cell = (lattice, positions, numbers)

    symmetry = spglib.get_symmetry(cell, symprec=tol)

    rotations, translations = symmetry['rotations'], symmetry['translations']

    return rotations, translations


def matching_coordinates(position, comparators, tol=1e-8):
    """Get the indices of all points in a comparator list that are
    equal to a given position (with a tolerance), taking into
    account periodic boundary conditions (adaptation from Pymatgen).

    This will only accept a Cartesian coordinate scheme.
    TODO: merge this with matching_sites.

    Parameters
    ----------
    position : list (3,)
        Fractional coordinate to compare to list.
    comparators : list (3, n)
        Fractional coordinates to compare against.
    tol : float
        Absolute tolerance.

    Returns
    -------
    match : list (n,)
        Indices of matches.
    """
    if len(comparators) == 0:
        return []

    fdist = comparators - position[None, :]
    match = np.where((np.abs(fdist) < tol).all(axis=1))[0]

    return match


def get_unique_coordinates(atoms, axis=2, tag=False, tol=1e-3):
    """Return unique coordinate values of a given atoms object
    for a specified axis.

    Parameters
    ----------
    atoms : object
        Atoms object to search for unique values along.
    axis : int (0, 1, or 2)
        Look for unique values along the x, y, or z axis.
    tag : bool
        Assign ASE-like tags to each layer of the slab.
    tol : float
        The tolerance to search for unique values within.

    Returns
    -------
    values : ndarray (n,)
        Array of unique positions in fractional coordinates.
    """
    positions = (atoms.get_scaled_positions()[:, axis] + tol) % 1
    positions -= tol

    values = [positions[0]]
    for d in positions[1:]:
        if not np.isclose(d, values, atol=tol, rtol=tol).any():
            values += [d]
    values = np.sort(values)

    if tag:
        tags = []
        for p in positions:
            close = np.isclose(p, values[::-1], atol=tol, rtol=tol)
            tags += [np.where(close)[0][0] + 1]
        atoms.set_tags(tags)

    return values


def get_affine_operations(atoms):
    """Return the affine operations of a unit cell."""
    rotations, translations = get_symmetry(atoms)

    operations = np.zeros((rotations.shape[0], 4, 4))
    for i, rot in enumerate(rotations):
        operations[i, :3, :3] = rot
        operations[i, -1, :3] = translations[i]
        operations[i, -1, -1] = 1

    return operations
