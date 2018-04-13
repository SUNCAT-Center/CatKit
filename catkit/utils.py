from . import Gratoms
from scipy.spatial import Voronoi
from scipy.linalg import lstsq
from ase.data import covalent_radii as radii
import numpy as np
from numpy.linalg import norm
import spglib
import os
import contextlib
import networkx as nx
import networkx.algorithms.isomorphism as iso


def trilaterate(centers, r):
    """Find the intersection of two or three spheres. In the case
    of two sphere intersection, the z-coordinate is assumed to be
    an intersection of a plane whose normal is aligned with the
    points and perpendicular to the positive z coordinate.

    Parameters:
    -----------
    centers : list or ndarray (n,)
        Array of values to have a average taken from.
    r : int
        Number of values to take an average with.

    Returns:
    --------
    intersection : ndarray (3,)
        The point where all spheres/planes intersect.
    """
    if len(r) == 1:
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

    if len(r) == 3:
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

    Parameters:
    -----------
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

    Parameters:
    -----------
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

    Parameters:
    -----------
    atoms : atoms object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.
    r : float
        Radius of the spheres to expand around each atom.

    Returns:
    --------
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

    Parameters:
    -----------
    atoms : atoms object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.
    r : float
        Radius of the spheres to expand around each atom.

    Returns:
    --------
    connectivity : ndarray (n, n)
        Number of edges formed between atoms in a system.
    """
    index, coords, offsets = expand_cell(atoms, r)
    L = int(len(offsets) / 2)
    oind = np.arange(L * len(atoms), (L + 1) * len(atoms))

    voronoi = Voronoi(coords)
    points = voronoi.ridge_points

    connectivity = np.zeros((len(atoms), len(atoms)), dtype=int)
    for i, n in enumerate(oind):
        p = points[np.where(points == n)[0]]
        edges = np.sort(index[p])

        unique_edge, edge_counts = np.unique(
            edges,
            return_counts=True,
            axis=0,
        )

        for j, edge in enumerate(unique_edge):
            u, v = edge
            connectivity[u][v] += edge_counts[j] / 2
            connectivity[v][u] += edge_counts[j] / 2

    return connectivity


def get_cutoff_neighbors(atoms, cutoff, atol=1e-8):
    """Return the connectivity matrix from a simple radial cutoff.
    Multi-bonding occurs through periodic boundary conditions.

    Parameters:
    -----------
    atoms : atoms object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.
    cutoff : float
        Cutoff radius to use when determining neighbors.
    atol: float
        Absolute tolerance to use when computing distances.

    Returns:
    --------
    connectivity : ndarray (n, n)
        Number of edges formed between atoms in a system.
    """
    cutoff += atol

    index, coords = expand_cell(atoms, cutoff * 2.0)[:2]

    connectivity = np.zeros((len(atoms), len(atoms)))
    for u, center in enumerate(atoms.positions):
        norm = np.linalg.norm(coords - center, axis=1)
        neighbors = index[np.where(norm < cutoff)]

        unique, counts = np.unique(neighbors, return_counts=True)
        connectivity[u, unique] = counts
        # Remove self-interaction
        connectivity[u, u] -= 1

    return connectivity


def get_neighbors(atoms, points=None, cutoff_matrix=None):
    """Returns the neighboring atoms within a specified cutoff matrix
    criteria for requested points.

    Use of the cutoff matrix provides more fine-tuned control
    over the interaction parameters.

    Parameters:
    -----------
    atoms : atoms object
        Atoms object to return.
    points : list (N,) or None
        Points to locate neighboring points of. If not provided, all
        atom points in the atoms-object will be used.
    cutoff_matrix : ndarray (96, 96) or None
        A matrix of interaction distances for the first 96 elements
        of the periodic table. These interactions are separated into
        individual i, j interactions. If None, defaults from ASE
        will be used.

    Returns:
    --------
    neighbors : dict
        Keys of each point and an array of each neighboring atoms index.
    """
    if cutoff_matrix is None:

        r = radii.copy()
        # TODO: Develop an SVM to parameterize this for me
        # Will need reliable training data for an supervised approach
        metals = [
            [47, 1.1],  # Ag
            [79, 1.2],  # Au
            [29, 1.1],  # Cu
            [77, 1.1],  # Ir
            [46, 1.1],  # Pd
            [78, 1.2],  # Pt
            [45, 1.1],  # Rh
        ]

        adsorbates = [[1, 1.0], [6, 1.0], [7, 1.0], [8, 1.0], [16, 1.0]]

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

                    neighbors[a] = np.concatenate((neighbors[a], i))

    return neighbors


def get_primitive_cell(atoms, tol=1e-8):
    """Atoms object interface with spglib primitive cell finder:
    https://atztogo.github.io/spglib/python-spglib.html#python-spglib

    Parameters:
    -----------
    atoms : object
        Atoms object to search for a primitive unit cell.
    tol : float
        Tolerance for floating point rounding errors.

    Returns:
    --------
    primitive cell : object
        The primitive unit cell returned by spglib if one is found.
    """
    lattice = atoms.cell
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()

    cell = (lattice, positions, numbers)

    _lattice, _positions, _numbers = spglib.find_primitive(cell, symprec=tol)

    atoms = Gratoms(symbols=_numbers, cell=_lattice, pbc=atoms.pbc)
    atoms.set_scaled_positions(_positions)

    return atoms


def get_symmetry(atoms, tol=1e-8):
    """Atoms object interface with spglib symmetry finder:
    https://atztogo.github.io/spglib/python-spglib.html#python-spglib

    Parameters:
    -----------
    atoms : object
        Atoms object to search for symmetric structures of.
    tol : float
        Tolerance for floating point rounding errors.

    Returns:
    --------
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

    Parameters:
    -----------
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

    Returns:
    --------
    values : ndarray (n,)
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


def plane_normal(xyz):
    """Return the surface normal vector to a plane of best fit.

    Parameters:
    -----------
    xyz : ndarray (n, 3)
        3D points to fit plane to.

    Returns:
    --------
    vec : ndarray (1, 3)
        Unit vector normal to the plane of best fit.
    """
    A = np.c_[xyz[:, 0], xyz[:, 1], np.ones(xyz.shape[0])]
    vec, _, _, _ = lstsq(A, xyz[:, 2])
    vec[2] = -1.0

    vec /= -np.linalg.norm(vec)

    return vec


def connectivity_to_edges(connectivity):
    """Convert a Numpy connectivity matrix into a list of
    NetworkX compatible edges
    """
    edges = []
    for i, c in enumerate(connectivity):
        lower_diagonal = c[:i]
        for j, v in enumerate(lower_diagonal):
            edges += [(int(i), int(j), 1)] * int(v)

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
