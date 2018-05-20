from . import defaults
from catkit import Gratoms
from scipy.spatial import Voronoi
from scipy.linalg import lstsq
from ase.data import chemical_symbols as sym
import numpy as np
from numpy.linalg import norm, solve
import networkx.algorithms.isomorphism as iso
import networkx as nx
import contextlib
import spglib
import os
import re


def trilaterate(centers, r):
    """Find the intersection of two or three spheres. In the case
    of two sphere intersection, the z-coordinate is assumed to be
    an intersection of a plane whose normal is aligned with the
    points and perpendicular to the positive z coordinate.

    Parameters
    ----------
    centers : list or ndarray (n,)
        Array of values to have a average taken from.
    r : int
        Number of values to take an average with.

    Returns
    -------
    intersection : ndarray (3,)
        The point where all spheres/planes intersect.
    """
    if len(r) > 3:
        raise ValueError('Cannot trilaterate more than 3 points')
    elif len(r) == 1:
        return centers[0] + [0, 0, r[0]]

    vec1 = centers[1] - centers[0]
    uvec1 = vec1 / norm(vec1)
    d = norm(vec1)

    if len(r) == 2:
        x0 = (d**2 - r[0]**2 + r[1]**2)
        x = d - x0 / (2 * d)
        a = np.sqrt(4 * d**2 * r[1]**2 - x0**2)
        z = 0.5 * (1 / d) * a
        if np.isnan(z):
            z = 0.01
        h = [0, 0, z]
        intersection = centers[0] + uvec1 * x + h

    elif len(r) == 3:
        vec2 = centers[2] - centers[0]
        i = np.dot(uvec1, vec2)
        vec2 = vec2 - i * uvec1

        uvec2 = vec2 / norm(vec2)
        uvec3 = np.cross(uvec1, uvec2)
        j = np.dot(uvec2, vec2)

        x = (r[0]**2 - r[1]**2 + d**2) / (2 * d)
        y = (r[0]**2 - r[2]**2 - 2 * i * x + i**2 + j**2) / (2 * j)
        z = np.sqrt(r[0]**2 - x**2 - y**2)
        if np.isnan(z):
            z = 0.01
        intersection = centers[0] + x * uvec1 + y * uvec2 + z * uvec3

    return intersection


@contextlib.contextmanager
def cd(path):
    """Does path management: if the path doesn't exists, create it
    otherwise, move into it until the indentation is broken.

    Parameters
    ----------
    path : str
        Directory path to create and change into.
    """
    cwd = os.getcwd()
    try:
        if not os.path.exists(path):
            os.makedirs(path)
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)


def rmean(x, N=5):
    """Calculate the running mean of array x for N instances.

    Parameters
    ----------
    x : list or ndarray (n,)
        Array of values to have a average taken from.
    N : int
        Number of values to take an average with.

    Returns:
    --------
    rmean : ndarray (n + 1,)
        Mean value of the running average.
    """
    length = len(x)
    if length < N:
        N = length

    cumsum = np.cumsum(np.insert(x, 0, 0))
    mean = (cumsum[N:] - cumsum[:-N]) / float(N)

    return mean


def expand_cell(atoms, r=6):
    """Return Cartesian coordinates atoms within a supercell
    which contains spheres of specified cutoff radius around
    all atom positions.

    Parameters
    ----------
    atoms : atoms object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.
    r : float
        Radius of the spheres to expand around each atom.

    Returns
    -------
    index : ndarray of int
        Indices associated with the original unit cell positions.
    coords : ndarray of (3,) array
        Cartesian coordinates associated with positions in the
        super-cell.
    """
    cell = atoms.get_cell()
    recp_len = np.diag(np.linalg.pinv(cell))
    nmax = float(r) * recp_len + 0.01

    pbc = atoms.get_pbc()
    low = np.floor(-nmax * pbc)
    high = np.ceil(nmax * pbc + 1)

    offsets = np.mgrid[low[0]:high[0], low[1]:high[1], low[2]:high[2]].T
    ncell = np.prod(offsets.shape[:-1])

    cart = np.dot(offsets, atoms.cell)

    coords = atoms.positions[None, None, None, :, :] + cart[:, :, :, None, :]

    index = np.arange(len(atoms))[None, :].repeat(ncell, axis=0).flatten()
    coords = coords.reshape(np.prod(coords.shape[:-1]), 3)
    offsets = offsets.reshape(ncell, 3)

    return index, coords, offsets


def get_voronoi_neighbors(atoms, r=8):
    """Return the connectivity matrix from the Voronoi
    method. Multi-bonding occurs through periodic boundary conditions.

    Parameters
    ----------
    atoms : atoms object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.
    r : float
        Radius of the spheres to expand around each atom.

    Returns
    -------
    connectivity : ndarray (n, n)
        Number of edges formed between atoms in a system.
    """
    index, coords, offsets = expand_cell(atoms, r)
    L = int(len(offsets) / 2)
    oind = np.arange(L * len(atoms), (L + 1) * len(atoms))

    voronoi = Voronoi(coords)
    points = voronoi.ridge_points

    connectivity = np.zeros((len(atoms), len(atoms)))
    for i, n in enumerate(oind):
        p = points[np.where(points == n)[0]]
        edges = np.sort(index[p])

        unique_edge, edge_counts = np.unique(
            edges,
            return_counts=True,
            axis=0)

        for j, edge in enumerate(unique_edge):
            u, v = edge
            connectivity[u][v] += edge_counts[j]
            connectivity[v][u] += edge_counts[j]

    connectivity /= 2

    return connectivity.astype(int)


def get_cutoff_neighbors(atoms, cutoff=None, atol=1e-8):
    """Return the connectivity matrix from a simple radial cutoff.
    Multi-bonding occurs through periodic boundary conditions.

    Parameters
    ----------
    atoms : atoms object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.
    cutoff : float
        Cutoff radius to use when determining neighbors.
    atol: float
        Absolute tolerance to use when computing distances.

    Returns
    -------
    connectivity : ndarray (n, n)
        Number of edges formed between atoms in a system.
    """
    cov_radii = defaults.get('covalent_radii')
    numbers = atoms.numbers
    index, coords = expand_cell(atoms, 4)[:2]

    if cutoff is None:
        radii = cov_radii[numbers]
        radii = np.repeat(radii[:, None], len(index) / len(radii), axis=1)
        radii = radii.T.flatten()
    else:
        radii = np.ones_like(index) * cutoff

    connectivity = np.zeros((len(atoms), len(atoms)), dtype=int)
    for i, center in enumerate(atoms.positions):
        norm = np.linalg.norm(coords - center, axis=1)
        boundry_distance = norm - radii
        if cutoff is None:
            boundry_distance -= cov_radii[numbers[i]]
        neighbors = index[np.where(boundry_distance < atol)]

        unique, counts = np.unique(neighbors, return_counts=True)
        connectivity[i, unique] = counts
        # Remove self-interaction
        connectivity[i, i] -= 1

    return connectivity


def get_primitive_cell(atoms, tol=1e-8):
    """Atoms object interface with spglib primitive cell finder:
    https://atztogo.github.io/spglib/python-spglib.html#python-spglib

    Parameters
    ----------
    atoms : object
        Atoms object to search for a primitive unit cell.
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
    cell = spglib.find_primitive(cell, symprec=tol)

    if cell is None:
        return None

    _lattice, _positions, _numbers = cell
    atoms = Gratoms(symbols=_numbers, cell=_lattice, pbc=atoms.pbc)
    atoms.set_scaled_positions(_positions)

    return atoms


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

    return spglib.get_symmetry(cell, symprec=tol)


def get_unique_coordinates(atoms, axis=2, direct=True, tag=False, tol=1e-5):
    """Return unique coordinate values of a given atoms object
    for a specified axis.

    Parameters
    ----------
    atoms : object
        Atoms object to search for unique values along.
    axis : int
        Value of 0, 1, or 2 associated with x, y, and z coordinates.
    direct : bool
        Whether to use direct coordinates or Cartesian.
    tag : bool
        Assign ase-like tags to each layer of the slab.
    tol : float
        The tolerance to search for unique values within.

    Returns
    -------
    values : ndarray (n,)
        Array of unique values.
    """
    if direct:
        positions = (atoms.get_scaled_positions() + tol) % 1
        positions -= tol
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


def plane_normal(xyz):
    """Return the surface normal vector to a plane of best fit.

    Parameters
    ----------
    xyz : ndarray (n, 3)
        3D points to fit plane to.

    Returns
    -------
    vec : ndarray (1, 3)
        Unit vector normal to the plane of best fit.
    """
    A = np.c_[xyz[:, 0], xyz[:, 1], np.ones(xyz.shape[0])]
    vec, _, _, _ = lstsq(A, xyz[:, 2])
    vec[2] = -1.0

    vec /= -np.linalg.norm(vec)

    return vec


def connectivity_to_edges(connectivity, indices=None):
    """Convert a Numpy connectivity matrix into a list of
    NetworkX compatible edges
    """
    if indices is None:
        indices = np.arange(connectivity.shape[0], dtype=int)

    edges = []
    for i, c in enumerate(connectivity):
        lower_diagonal = c[:i]
        for j, v in enumerate(lower_diagonal):
            edges += [(indices[i], indices[j], 1)] * int(v)

    return edges


def isomorphic_molecules(graph0, graph1):
    """Check whether two molecule graphs are isomorphic."""
    em = iso.numerical_edge_match('bonds', 1)
    nm = iso.numerical_node_match('number', 1)

    isomorphic = nx.is_isomorphic(graph0, graph1, edge_match=em, node_match=nm)

    return isomorphic


def to_gratoms(atoms):
    """Convert and atom object to a gratoms object."""

    gratoms = Gratoms(
        numbers=atoms.numbers,
        positions=atoms.positions,
        pbc=atoms.pbc,
        cell=atoms.cell)

    if atoms.constraints:
        gratoms.set_constraint(atoms.constraints)

    return gratoms


def get_atomic_numbers(formula, return_count=False):
    """Return the atomic numbers associated with a chemical formula.

    Parameters
    ----------
    formula : string
        A chemical formula to parse into atomic numbers.
    return_count : bool
        Return the count of each element in the formula.

    Returns
    -------
    numbers : ndarray (n,)
        Element numbers in associated species.
    counts : ndarray (n,)
        Count of each element in a species.
    """
    parse = re.findall('[A-Z][a-z]?|[0-9]+', formula)

    values = {}
    for i, e in enumerate(parse):
        if e.isdigit():
            values[parse[i - 1]] += int(e) - 1
        else:
            if e not in values:
                values[e] = 1
            else:
                values[e] += 1

    numbers = np.array([sym.index(k) for k in values.keys()])
    srt = np.argsort(numbers)
    numbers = numbers[srt]

    if return_count:
        counts = np.array([v for v in values.values()])[srt]

        return numbers, counts

    return numbers


def get_reference_energies(species, energies):
    """Get reference energies for the elements in a set of molecules.

    Parameters
    ----------
    species : list (n,)
        Chemical formulas for each molecular species.
    energies : list (n,)
        Total energies associated with each species.

    Returns
    -------
    elements : ndarray (n,)
        Atomic elements associated with all species.
    references : ndarray (n,)
        Reference energies associated with each element.
    """
    if not isinstance(energies, np.ndarray):
        energies = np.array(energies)

    A = np.zeros((len(species), len(species)))
    elements = np.zeros(len(species), dtype=int)
    n = 0

    # Construct the elements array as they appear
    for i, s in enumerate(species):
        num, cnt = get_atomic_numbers(s, True)

        for j in num[~np.in1d(num, elements)]:
            elements[n] = j
            n += 1

        A[i][np.in1d(elements, num)] = cnt

    references = solve(A, energies)
    srt = np.argsort(elements)
    references = references[srt]
    elements = elements[srt]

    return elements, references


def get_unique_xy(xyz_coords, cutoff=0.1):
    """Return the unique coordinates of an atoms object
    for the requrested atoms indices. Z-coordinates are projected
    to maximum z-coordinate by default.

    Parameters
    ----------
    xyz_coords : ndarray (n, 3)
        Cartesian coordinates to identify unique xy positions from.
    cutoff : float
        Distance in Angstrons to consider xy-coordinate unique within.

    Returns
    -------
    xy_pos : ndarray (m, 3)
        Unique xy coordinates projected onto a maximal z coordinate.
    """
    xyz_coords[:, -1] = np.max(xyz_coords[:, -1])

    xy_copies = []
    for i, p in enumerate(xyz_coords[:, :-1]):
        if i in xy_copies:
            continue

        dis = xyz_coords[:, :-1][:, None] - p
        match = np.where(abs(dis).sum(axis=2).T < cutoff)[1]
        xy_copies += match[1:].tolist()

    xyz_coords = np.delete(xyz_coords, xy_copies, axis=0)

    return xyz_coords


def parse_slice(slice_name):
    """Return a correctly parsed slice from input of
    varying types.
    """
    if isinstance(slice_name, (slice)):
        _slice = slice_name

    elif isinstance(slice_name, type(None)):
        _slice = slice(None)

    elif isinstance(slice_name, int):
        i = int(slice_name)
        _slice = slice(i, i + 1)

    elif isinstance(slice_name, str):
        if slice_name.isdigit():
            i = int(slice_name)
            _slice = slice(i, i + 1)

        else:
            split = slice_name.split(':')
            split = [int(_) if _.lstrip('-').isdigit()
                     else None for _ in split]
            _slice = slice(*split)

    return _slice
