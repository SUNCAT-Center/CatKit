from scipy.spatial import Voronoi
from ase.data import covalent_radii as radii
import numpy as np
import spglib
from ase import Atoms


def expand_cell(atoms, r=6):
    """ Return Cartesian coordinates atoms within a supercell
    which contains spheres of specified cutoff radius around
    all atom positions.

    Args:
      atoms: ASE atoms-object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.

      r: float
        Radius of the spheres to expand around each atom.

    Returns:
      index: ndarray of int
        Indices associated with the original unit cell positions.

      coords: ndarray of (3,) array
        Cartesian coordinates associated with positions in the
        supercell.
    """

    cell = atoms.get_cell()
    recp_len = np.diag(np.linalg.pinv(cell))
    nmax = float(r) * recp_len + 0.01

    pbc = atoms.get_pbc()
    low = np.floor(-nmax * pbc)
    high = np.ceil(nmax * pbc + 1)

    arange = np.arange(low[0], high[0])
    brange = np.arange(low[1], high[1])
    crange = np.arange(low[2], high[2])

    arange = arange[:, None] * np.array([1, 0, 0])[None, :]
    brange = brange[:, None] * np.array([0, 1, 0])[None, :]
    crange = crange[:, None] * np.array([0, 0, 1])[None, :]

    images = arange[:, None, None] + \
             brange[None, :, None] + \
             crange[None, None, :]

    cart_images = np.dot(images, cell)

    coords = atoms.positions[:, None, None, None, :] + \
             cart_images[None, :, :, :, :]

    index = np.where(coords.sum(axis=4))[0]
    index = np.append([0], index)

    coords = coords.flatten()
    n = int(coords.shape[0] / 3)
    coords = coords.reshape((n, 3))

    return index, coords


def get_voronoi_neighbors(atoms, r=10):
    """ Return the nearest-neighbors list from the Voronoi
    method.

    Args:
      atoms: ASE atoms-object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.

      r: float
        Radius of the spheres to expand around each atom.

    Returns:
      neighbors: dict
        Edge tuples notation denoting bonds between the atoms of
        the corresponding indices and values of the number of
        those bonds.

        Multi-bonding occurs through periodic boundary conditions.
    """

    index, coords = expand_cell(atoms, r)

    oind = np.empty(len(atoms))
    for i, A in enumerate(atoms.positions):
        for j, B in enumerate(coords):
            if np.allclose(A, B):
                oind[i] = j
                break

    voronoi = Voronoi(coords)
    edges = {}
    indices = np.zeros(len(atoms))
    for nn, vind in voronoi.ridge_dict.items():
        if -1 in vind:
            continue

        for n in nn:
            if n in oind:
                indices[index[n]] += 1
                e = tuple(sorted((index[nn[0]], index[nn[1]])))
                if e not in edges:
                    edges[e] = 0.5
                else:
                    edges[e] += 0.5

    return indices, edges


def get_cutoff_neighbors(atoms, cutoff=4, atol=1e-8):
    """ I can haz documentation?
    """

    cutoff = cutoff + atol

    # TODO fix bug with 110
    index, coords = expand_cell(atoms, cutoff * 2.0)

    indices = np.zeros(len(atoms))
    edges = {}
    for i, center in enumerate(atoms.positions):
        norm = coords - center
        d2 = (norm ** 2).sum(axis=1)
        srt = np.argsort(d2)

        d2 = d2[srt][1:]
        nindex = index[srt][1:]

        inner = np.where(d2 <= cutoff ** 2)[0]

        for j in nindex[inner]:
            indices[i] += 1
            e = tuple(sorted((i, j)))
            if e not in edges:
                edges[e] = 0.5
            else:
                edges[e] += 0.5

    return indices, edges


def get_neighbors(
        atoms,
        points=None,
        cutoff_matrix=None):
    """ Returns the neighboring atoms within a specified cutoff matrix
    criteria for requested points.

    Use of the cutoff matrix provides more fine-tuned control
    over the interaction parameters.

    Args:
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
            [45, 1.1],  # Rh
        ]

        adsorbates = [
            [1, 1.0],
            [6, 1.0],
            [7, 1.0],
            [8, 1.0],
            [16, 1.0]
        ]

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

    Args:
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

    Args:
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
        tag=False,
        tol=1e-5):
    """ Return unique coordinate values of a given atoms object
    for a specified axis.

    Args:
      atoms: ASE atoms-object
        Atoms object to search for unique values along.

      axis: int
        Value of 0, 1, or 2 associated with x, y, and z coordinates.

      direct: bool
        Whether to use direct coordinates or Cartesian.

      tag: bool
        Assign ase-like tags to each layer of the slab.

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
    values = np.sort(values)

    if tag:
        tags = []
        for p in positions[:, axis]:
            close = np.isclose(p, values[::-1], rtol=tol)
            tags += [np.where(close)[0][0] + 1]
        atoms.set_tags(tags)

    return values
