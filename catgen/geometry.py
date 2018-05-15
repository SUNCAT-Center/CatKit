import numpy as np


def matching_sites(position, comparators, tol=1e-8):
    """Get the indices of all points in a comparator list that are
    equal to a given position (with a tolerance), taking into
    account periodic boundary conditions (adaptation from Pymatgen).

    Parameters:
    -----------
    position : list (3,)
        Fractional coordinate to compare to list.
    comparators : list (3, n)
        Fractional coordinates to compare against.
    tol : float
        Absolute tolerance.

    Returns:
    --------
    match : list (n,)
        Indices of matches.
    """
    if len(comparators) == 0:
        return []

    fdist = comparators - position
    fdist -= np.round(fdist)
    match = np.where((np.abs(fdist) < tol).all(axis=1))[0]

    return match


# def build_xyz(self):
#     """ Build xyz representation from z-matrix"""

#     coords = self.atomcoords[-1]
#     self.newcoords = np.zeros((len(coords), 3))

#     for i in range(len(coords)):
#         self.newcoords[i] = self._calc_position(i)

#     self.atomcoords[-1] = self.newcoords


def _calc_position(self, i, connectivity, angleconnectivity, dihedralconnectivity):
    """Calculate position of another atom based on internal coordinates"""

    if i > 1:
        j = connectivity[i]
        k = angleconnectivity[i]
        l = dihedralconnectivity[i]

        # # Prevent doubles
        # if k == l and i > 0:
        #     for idx in range(1, len(self.connectivity[:i])):
        #         if self.connectivity[idx] in [i, j, k] and not idx in [i, j, k]:
        #             l = idx
        #             break

        avec = self.newcoords[j]
        bvec = self.newcoords[k]

        dst = self.distances[i]
        ang = self.angles[i] * np.pi / 180.0

        if i == 2:
            # Third atom will be in same plane as first two
            tor = 90.0 * np.pi / 180.0
            cvec = np.array([0, 1, 0])
        else:
            # Fourth + atoms require dihedral (torsional) angle
            tor = self.dihedrals[i] * np.pi / 180.0
            cvec = self.newcoords[l]

        v1 = avec - bvec
        v2 = avec - cvec

        n = np.cross(v1, v2)
        nn = np.cross(v1, n)

        n /= np.norm(n)
        nn /= np.norm(nn)

        n *= -np.sin(tor)
        nn *= np.cos(tor)

        v3 = n + nn
        v3 /= np.norm(v3)
        v3 *= dst * np.sin(ang)

        v1 /= np.norm(v1)
        v1 *= dst * np.cos(ang)

        position = avec + v3 - v1

    elif i == 1:
        # Second atom dst away from origin along Z-axis
        j = self.connectivity[i]
        dst = self.distances[i]
        position = np.array(
            [self.newcoords[j][0] + dst,
             self.newcoords[j][1],
             self.newcoords[j][2]])

    elif i == 0:
        # First atom at the origin
        position = np.array([0, 0, 0])

    return position
