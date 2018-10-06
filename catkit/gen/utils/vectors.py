import numpy as np
import scipy


def get_reciprocal_vectors(atoms):
    """Return the reciprocal lattice vectors to a atoms unit cell."""
    rotation1 = np.roll(atoms.cell, -1, axis=0)
    rotation2 = np.roll(atoms.cell, 1, axis=0)
    normal_vectors = np.cross(rotation1, rotation2)

    return normal_vectors


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
    vec, _, _, _ = scipy.linalg.lstsq(A, xyz[:, 2])
    vec[2] = -1.0

    vec /= -np.linalg.norm(vec)

    return vec


def get_basis_vectors(coordinates):
    """Return a set of basis vectors for a given array of
    3D coordinates.

    Parameters
    ----------
    coordinates : array_like (3, 3) | (2, 3)
        Cartesian coordinates to determine the basis of. If
        only 2 positions are given 3rd is chosen as the positive
        y-axis.

    Returns
    -------
    basis_vectors : ndarray (3, 3)
        Automatically generated basis vectors from the given
        positions.
    """
    if len(coordinates) == 3:
        c0, c1, c2 = coordinates
    else:
        c0, c1 = coordinates
        c2 = np.array([0, 1, 0])

    basis1 = c0 - c1
    basis2 = np.cross(basis1, c0 - c2)
    basis3 = np.cross(basis1, basis2)

    basis_vectors = np.vstack([basis1, basis2, basis3])
    basis_vectors /= np.linalg.norm(
        basis_vectors, axis=1, keepdims=True)

    return basis_vectors
