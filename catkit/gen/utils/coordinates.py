import numpy as np


def expand_cell(
        coordinates,
        cell,
        pbc=True,
        centers=None,
        cutoff=6.0):
    """Return Cartesian coordinates atoms within a supercell
    which contains spheres of a given cutoff radius. One sphere
    for each atom inside of the given unit cell.

    Parameters
    ----------
    coordinates : ndarray (N, 3)
        Cartesian coordinates of atomic positions.
    cell : ndarray (3, 3)
        Unit cell dimensions
    pbc : bool | array_like (3,)
        Periodic boundary conditions.
    cutoff : float
        Radius of the spheres to use for expanding the unit cell.

    Returns
    -------
    indices : ndarray (M,)
        Indices associated with the original unit cell atom index.
    positions : ndarray (M, 3)
        Cartesian coordinates associated with positions in the
        supercell.
    offsets : ndarray (M, 3)
        Integer offsets of each unit cell.
    """
    pbc = np.asarray(pbc, dtype=bool)
    if centers is None:
        centers = np.linalg.solve(cell.T, coordinates.T).T

    basis_lengths = np.linalg.norm(cell, axis=1)
    reciprocal_cell = np.linalg.inv(cell).T
    reciprocal_lengths = np.linalg.norm(reciprocal_cell, axis=1)

    lengths = pbc * cutoff * reciprocal_lengths + 0.01
    min_pad = np.floor(centers - lengths).min(axis=0)
    max_pad = np.ceil(centers + lengths).max(axis=0)
    padding = np.array([min_pad.astype(int), max_pad.astype(int)]).T

    offsets = np.mgrid[[slice(*_) for _ in padding]].T

    ncell = np.prod(offsets.shape[:-1])
    indices = np.tile(np.arange(len(coordinates)), ncell)

    cartesian_cells = np.dot(offsets, cell)
    positions = cartesian_cells[:, :, :, None, :] + \
        coordinates[None, None, None, :, :]
    positions = positions.reshape(-1, 3)

    offsets = np.tile(offsets, len(coordinates)).reshape(-1, 3)

    return indices, positions, offsets


def trilaterate(centers, r, zvector=None):
    """Find the intersection of two or three spheres. In the case
    of two sphere intersection, the z-coordinate is assumed to be
    an intersection of a plane whose normal is aligned with the
    points and perpendicular to the positive z-coordinate.

    If more than three spheres are supplied, the centroid of the
    points is returned (no radii used).

    Parameters
    ----------
    centers : list | ndarray (n, 3)
        Cartesian coordinates representing the center of each sphere
    r : list | ndarray (n,)
        The radii of the spheres.
    zvector : ndarray (3,)
        The vector associated with the upward direction for under-specified
        coordinations (1 and 2).

    Returns
    -------
    intersection : ndarray (3,)
        The point where all spheres/planes intersect.
    """
    if zvector is None:
        zvector = np.array([0, 0, 1])

    if len(r) == 1:
        return centers[0] + r[0] * zvector
    elif len(r) > 3:
        centroid = np.sum(centers, axis=0) / len(centers)
        centroid += np.mean(r) / 2 * zvector
        return centroid

    vec1 = centers[1] - centers[0]
    uvec1 = vec1 / np.linalg.norm(vec1)
    d = np.linalg.norm(vec1)

    if len(r) == 2:
        x0 = (d**2 - r[0]**2 + r[1]**2)
        x = d - x0 / (2 * d)
        a = np.sqrt(4 * d**2 * r[1]**2 - x0**2)
        z = 0.5 * (1 / d) * a
        if np.isnan(z):
            z = 0.01
        h = z * zvector
        intersection = centers[0] + uvec1 * x + h

    elif len(r) == 3:
        vec2 = centers[2] - centers[0]
        i = np.dot(uvec1, vec2)
        vec2 = vec2 - i * uvec1

        uvec2 = vec2 / np.linalg.norm(vec2)
        uvec3 = np.cross(uvec1, uvec2)
        j = np.dot(uvec2, vec2)

        x = (r[0]**2 - r[1]**2 + d**2) / (2 * d)
        y = (r[0]**2 - r[2]**2 - 2 * i * x + i**2 + j**2) / (2 * j)
        z = np.sqrt(r[0]**2 - x**2 - y**2)
        if np.isnan(z):
            z = 0.01
        intersection = centers[0] + x * uvec1 + y * uvec2 + z * uvec3

    return intersection


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


def matching_sites(position, comparators, tol=1e-8):
    """Get the indices of all points in a comparator list that are
    equal to a given position (with a tolerance), taking into
    account periodic boundary conditions (adaptation from Pymatgen).

    This will only accept a fractional coordinate scheme.

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

    fdist = comparators - position
    fdist -= np.round(fdist)
    match = np.where((np.abs(fdist) < tol).all(axis=1))[0]

    return match


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
    comparators : list (3, N)
        Fractional coordinates to compare against.
    tol : float
        Absolute tolerance.

    Returns
    -------
    match : list (N,)
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


def get_integer_enumeration(N=3, span=[0, 2]):
    """Return the enumerated array of a span of integer values.
    These enumerations are limited to the length N.

    For the default span of [0, 2], the enumeration equates to
    the corners of an N-dimensional hypercube.

    Parameters
    ----------
    N : int
        Length of enumerated lists to produce.
    span : list | slice
        The range of integers to be considered for enumeration.

    Returns
    -------
    enumeration : ndarray (M, N)
        Enumeration of the requested integers.
    """
    if not isinstance(span, slice):
        span = slice(*span)
    enumeration = np.mgrid[[span] * N].reshape(N, -1).T

    return enumeration
