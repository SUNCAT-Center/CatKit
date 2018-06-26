import numpy as np


def expand_cell(atoms, r=6):
    """Return Cartesian coordinates atoms within a supercell
    which contains spheres of specified cutoff radius around
    all atom positions.

    Parameters
    ----------
    atoms : Atoms object
        Atoms with the periodic boundary conditions and unit cell
        information to use.
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
    low = np.floor(-np.abs(nmax) * pbc)
    high = np.ceil(np.abs(nmax) * pbc + 1)

    offsets = np.mgrid[low[0]:high[0], low[1]:high[1], low[2]:high[2]].T

    ncell = np.prod(offsets.shape[:-1])
    cart = np.dot(offsets, atoms.cell)
    coords = atoms.positions[None, None, None, :, :] + cart[:, :, :, None, :]

    index = np.arange(len(atoms))[None, :].repeat(ncell, axis=0).flatten()
    coords = coords.reshape(np.prod(coords.shape[:-1]), 3)
    offsets = offsets.reshape(ncell, 3)

    return index, coords, offsets


def trilaterate(centers, r):
    """Find the intersection of two or three spheres. In the case
    of two sphere intersection, the z-coordinate is assumed to be
    an intersection of a plane whose normal is aligned with the
    points and perpendicular to the positive z-coordinate.

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
    uvec1 = vec1 / np.linalg.norm(vec1)
    d = np.linalg.norm(vec1)

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
