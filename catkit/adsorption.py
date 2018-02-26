from . import utils
from scipy.spatial import Delaunay
from scipy.linalg import circulant
from itertools import product
from numpy.linalg import pinv, norm
import numpy as np

# TODO:
# get_reduced_sites():
# Should return indices of the reduced sites
# 

class AdsorptionSites():
    """Adsorption site object."""

    def __init__(self, slab, r=8, mask=0.5, tol=1e-5):
        """Create an extended unit cell of the surface sites for
        use in identifying other sites.

        TODO: Determine if this can be made more efficient by
        removing the 'sites' dictionary.

        Parameters:
        -----------
        slab : atoms object
            The atoms object to manage adsorption sites for. Must contain
            surface atom identification.
        r : float
            Minimum basis vector length in Angstroms for creating extended
            unit cell.
        mask : float
            Percentage of fractional coordinates to ignore.
        tol : float
            Absolute tolerance for floating point errors.
        """
        index, coords, offsets = utils.expand_cell(slab, r)
        surface_atoms = slab.get_surface_atoms()
        if surface_atoms is None:
            raise ValueError('Slab must contain surface atoms')

        extended_top = np.where(np.in1d(index, surface_atoms))[0]

        self.tol = tol
        self.coordinates = coords[extended_top].tolist()
        self.connectivity = np.ones(extended_top.shape[0]).tolist()
        self.r1_topology = [[et] for et in extended_top]
        self.unique_index = []

        sites = self.get_higher_coordination_sites(
            coords[extended_top]
        )
        self.r2_topology = sites['top'][2]

        # Put data into array format
        for i, items in enumerate(sites.items()):
            k, v = items
            if k is 'top':
                continue
            coords, r1top, r2top = v

            self.connectivity += (np.ones(len(coords)) * (i + 1)).tolist()
            self.coordinates += coords
            self.r1_topology += r1top
            self.r2_topology += r2top

        self.coordinates = np.array(self.coordinates)
        self.connectivity = np.array(self.connectivity, dtype=int)
        self.r1_topology = np.array(self.r1_topology)
        self.r2_topology = np.array(self.r2_topology)
        self.frac_coords = np.dot(self.coordinates, pinv(slab.cell))
        self.offsets = offsets
        self.slab_index = index
        self.slab = slab

        xlim = offsets.T[0].max()
        ylim = offsets.T[1].max()

        screen = (
            self.frac_coords[:, 0] < xlim + 1 - (mask * xlim)) & \
            (self.frac_coords[:, 0] > -xlim + (mask * xlim)) & \
            (self.frac_coords[:, 1] < ylim + 1 - (mask * ylim)) & \
            (self.frac_coords[:, 1] > -ylim + (mask * ylim))

        self.mask = screen

    def get_higher_coordination_sites(
            self,
            top_coordinates,
    ):
        """Find all bridge and hollow sites (3-fold and 4-fold) given an
        input slab based Delaunay triangulation of surface atoms of a
        super-cell.

        Parameters:
        -----------
        top_coordinates : ndarray (n, 3)
            Cartesian coordinates for the top atoms of the unit cell.

        Returns:
        --------
        sites : dict of 3 lists
            Dictionary sites containing positions, points, and neighbor lists.
        """
        sites = {
            'top': [
                top_coordinates,
                [],
                [[] for _ in top_coordinates]
            ],
            'bridge': [[], [], []],
            'hollow': [[], [], []],
            '4fold': [[], [], []],
        }

        dt = Delaunay(sites['top'][0][:, :2])
        neighbors = dt.neighbors
        simplices = dt.simplices

        for i, corners in enumerate(simplices):
            cir = circulant(corners)
            edges = cir[:, 1:]

            # Inner angle of each triangle corner
            vec = sites['top'][0][edges.T] - sites['top'][0][corners]
            uvec = vec.T / norm(vec, axis=2).T
            angles = np.sum(uvec.T[0] * uvec.T[1], axis=1)

            # Angle types
            right = np.isclose(angles, 0)
            obtuse = (angles < -self.tol)

            rh_corner = corners[right]
            edge_neighbors = neighbors[i]

            if obtuse.any():
                # Assumption: All simplices with obtuse angles
                # are irrelevant boundaries.
                continue

            bridge = np.sum(sites['top'][0][edges], axis=1) / 2.0

            # Looping through corners allows for elimination of
            # redundant points, identification of 4-fold hollows,
            # and collection of bridge neighbors.
            for j, c in enumerate(corners):
                edge = sorted(edges[j])

                if edge in sites['bridge'][1]:
                    continue

                # Get the bridge neighbors (for adsorption vector)
                neighbor_simplex = simplices[edge_neighbors[j]]
                oc = list(set(neighbor_simplex) - set(edge))[0]

                # Right angles potentially indicate 4-fold hollow
                potential_hollow = edge + sorted([c, oc])
                if c in rh_corner:

                    if potential_hollow in sites['4fold'][1]:
                        continue

                    # Assumption: If not 4-fold, this suggests
                    # no hollow OR bridge site is present.
                    ovec = sites['top'][0][edge] - sites['top'][0][oc]
                    ouvec = ovec / norm(ovec)
                    oangle = np.dot(*ouvec)
                    oright = np.isclose(oangle, 0)
                    if oright:
                        sites['4fold'][0] += [bridge[j]]
                        sites['4fold'][1] += [potential_hollow]
                        sites['top'][2][c] += [oc]
                else:
                    sites['bridge'][0] += [bridge[j]]
                    sites['bridge'][1] += [edge]
                    sites['bridge'][2] += [[c, oc]]

                sites['top'][2][edge[0]] += [edge[1]]
                sites['top'][2][edge[1]] += [edge[0]]

            if not right.any():
                hollow = np.average(sites['top'][0][corners], axis=0)
                sites['hollow'][0] += [hollow]
                sites['hollow'][1] += [corners.tolist()]

        # For collecting missed bridge neighbors
        for s in sites['4fold'][1]:

            for edge in product(s[:2], s[2:]):
                edge = sorted(edge)
                i = sites['bridge'][1].index(edge)
                n, m = sites['bridge'][1][i], sites['bridge'][2][i]
                nn = (set(s) - set(n + m))
                sites['bridge'][2][i] += [list(nn)[0]]

        return sites

    def get_periodic_sites(self):
        """Return an index of the coordinates which are unique by
        periodic boundary conditions.

        Returns:
        --------
        periodic_match : ndarray (n,)
            Indices of the coordinates which are identical by
            periodic boundary conditions.
        """
        periodic_match = np.arange(self.frac_coords.shape[0])
        for i, j in enumerate(periodic_match):
            if i != j:
                continue

            new_match = matching_sites(self.frac_coords[i], self.frac_coords)
            periodic_match[new_match] = i

        return periodic_match

    def get_symmetric_sites(self):
        """Determine the symmetrically unique adsorption sites
        from a list of fractional coordinates.

        Returns:
        --------
        sites : dict of lists
            Dictionary of sites containing index of site
        """
        symmetry = utils.get_symmetry(self.slab, tol=self.tol)
        rotations = np.swapaxes(symmetry['rotations'], 1, 2)
        translations = symmetry['translations']
        affine = np.append(rotations, translations[:, None], axis=1)

        affine_points = np.insert(self.frac_coords, 3, 1, axis=1)
        operations = np.dot(affine_points, affine)

        symmetry_match = np.arange(self.frac_coords.shape[0])
        for i, j in enumerate(symmetry_match):
            if i != j:
                continue

            d = operations[i, :, None] - self.frac_coords
            d -= np.round(d)
            dind = np.where((np.abs(d) < self.tol).all(axis=2))[-1]
            symmetry_match[np.unique(dind)] = i

        return symmetry_match


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


def _get_adsorption_vectors(vslab, sites):
    """Returns the vectors representing the furthest distance from
    the neighboring atoms.

    (TODO: This input is complex and confusing. Would be nice to simplify.)

    Parameters:
    -----------
    vslab : object
        The virtual surface produced from find_adsorption_sites.
    sites : dict of 3 lists
        Dictionary of top, bridge, hollow, and 4-fold sites containing
        positions, points, and neighbor lists.

    Returns:
    --------
    sites : dict of 3 lists
        Dictionary of top, bridge, hollow, and 4-fold sites containing
        positions, points, and adsorption vector lists.
    """
    pos = vslab.positions
    for k, v in sites.items():
        coordinates, points, neighbors = v

        vectors = []
        for i, s in enumerate(coordinates):
            if len(neighbors):
                direct = pos[np.append(points[i], neighbors[i]).astype(int)]
            else:
                direct = pos[points[i]]

            vectors += [utils.plane_normal(direct)]

        sites[k][2] = np.array(vectors)

    return sites


def get_adsorption_sites(
        slab,
        periodic_reduced=True,
        symmetry_reduced=True,
        vectors=False,
        tol=1e-5,
):
    """Get the adsorption sites of a slab as defined by surface
    atom symmetries.

    Parameters:
    -----------
    slab : atoms object
        The slab to find adsorption sites for.
    symmetry_reduced : int
        Return the symmetrically unique sites only.
    vectors : bool
        Compute the adsorption vectors.

    Returns:
    --------
    sites : dict of 3 lists
        Dictionary of top, bridge, hollow, and 4-fold sites containing
        positions, points, and neighbor lists. If adsorption vectors
        are requested, the third list is replaced.
    """
    sites = find_adsorption_sites(
        slab=slab,
        trim=0.5,
        tol=tol,
    )

    if periodic_reduced:
        sites = get_reduced_sites(
            sites=sites,
            slab=slab,
            tol=tol,
        )

    if symmetry_reduced:
        sites = get_symmetric_sites(
            sites=sites,
            slab=slab,
            tol=tol,
        )

    if vectors:
        sites = _get_adsorption_vectors(sites)

    return sites
