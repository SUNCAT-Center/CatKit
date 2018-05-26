import numpy as np
from numpy.linalg import norm
from . import defaults


def matching_sites(position, comparators, tol=1e-8):
    """Get the indices of all points in a comparator list that are
    equal to a given position (with a tolerance), taking into
    account periodic boundary conditions (adaptation from Pymatgen).

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


def _get_basis_vectors(coordinates):
    if len(coordinates) == 3:
        c0, c1, c2 = coordinates
    else:
        c0, c1 = coordinates
        c2 = np.array([0, 1, 0])

    basis1 = c0 - c1
    basis2 = np.cross(basis1, c0 - c2)
    basis3 = np.cross(basis1, basis2)

    basis1 /= norm(basis1)
    basis2 /= norm(basis2)
    basis3 /= norm(basis3)

    basis_vectors = np.vstack([basis1, basis2, basis3])

    return basis_vectors


def _get_position(coordinates, basis=None,
                  distance=1, angle=109.47, dihedral=0):
    if basis is None:
        basis = _get_basis_vectors(coordinates)

    ang = np.deg2rad(angle)
    tor = np.deg2rad(dihedral)

    basis[1] *= -np.sin(tor)
    basis[2] *= np.cos(tor)

    vector = basis[1] + basis[2]
    vector /= norm(vector)

    vector *= distance * np.sin(ang)
    basis[0] *= distance * np.cos(ang)

    position = coordinates[0] + vector - basis[0]

    return position


def _branch_molecule(atoms, branch, base_root=None, adsorption=False):
    root, nodes = branch

    if len(nodes) == 0:
        return

    radii = defaults.get('covalent_radii')
    num = atoms.numbers[[root] + nodes]
    d = radii[num[0]] + radii[num[1:]]
    c0 = atoms[root].position

    angle = 120
    angle_mod = 0
    dihedral = np.array([0, 120, -120])

    if base_root is None:
        c1 = np.array([0, 0, -d[0]])

        if not adsorption:
            # Remove the first element for proper angle treatment
            m = nodes.pop(0)
            atoms[m].position = c1
            d = d[1:]
        else:
            # Move adsorption structures away from surface
            angle_mod = 25
    else:
        c1 = atoms[base_root].position

    n = len(nodes)
    if n == 1:
        if base_root is None and not adsorption:
            dihedral[0] = 120
        else:
            angle = 180
    elif n == 2:
        dihedral[1] = 180
    else:
        angle = 109.47

    for k in range(n):
        c = _get_position(
            [c0, c1],
            distance=d[k],
            angle=angle + angle_mod,
            dihedral=dihedral[k])
        atoms[nodes[k]].position = c

    return root
