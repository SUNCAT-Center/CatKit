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
