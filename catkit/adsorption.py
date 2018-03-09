from . import utils
from scipy.spatial import Delaunay
from scipy.linalg import circulant
from itertools import product
from numpy.linalg import pinv, norm
import numpy as np
from scipy.spatial import Voronoi
import matplotlib.pyplot as plt


class AdsorptionSites():
    """Adsorption site object."""

    def __init__(self, slab, r=8, tol=1e-5):
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
        self.r1_topology = [[i] for i in np.arange(len(extended_top))]
        self.index = index[extended_top]

        sites = self.get_higher_coordination_sites(coords[extended_top])
        self.r2_topology = sites['top'][2]

        # Put data into array format
        for i, items in enumerate(sites.items()):
            k, v = items
            if k is 'top':
                continue
            coords, r1top, r2top = v

            if k in ['hollow', '4fold']:
                r2top = [[] for _ in coords]

            self.connectivity += (np.ones(len(coords)) * (i + 1)).tolist()
            self.coordinates += coords
            self.r1_topology += r1top
            self.r2_topology += r2top

        self.coordinates = np.array(self.coordinates)
        self.connectivity = np.array(self.connectivity, dtype=int)
        self.r1_topology = np.array(self.r1_topology)
        self.r2_topology = np.array(self.r2_topology)
        self.frac_coords = np.dot(self.coordinates, pinv(slab.cell))
        self.slab = slab

        screen = (self.frac_coords[:, 0] > 0 - self.tol) & \
                 (self.frac_coords[:, 0] < 1 - self.tol) & \
                 (self.frac_coords[:, 1] > 0 - self.tol) & \
                 (self.frac_coords[:, 1] < 1 - self.tol)

        self.screen = screen

    def get_connectivity(self, screen=True):
        """Return the number of connections associated with each site."""
        if screen:
            return self.connectivity[self.screen]
        else:
            return self.connectivity

    def get_coordinates(self, screen=True):
        """Return the 3D coordinates associated with each site."""
        if screen:
            return self.coordinates[self.screen]
        else:
            return self.coordinates

    def get_topology(self, screen=True):
        """Return the 3D coordinates associated with each site."""
        topology = [self.index[top] for top in self.r1_topology]
        topology = np.array(topology)
        if screen:
            return topology[self.screen]
        else:
            return topology

    def get_higher_coordination_sites(
            self,
            top_coordinates,
            allow_obtuse=True):
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
            'top': [top_coordinates, [], [[] for _ in top_coordinates]],
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

            if obtuse.any() and not allow_obtuse:
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

            if not right.any() and not obtuse.any():
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

    def get_periodic_sites(self, screen=True):
        """Return an index of the coordinates which are unique by
        periodic boundary conditions.

        Parameters:
        -----------
        screen : bool
            Return only sites inside the unit cell.

        Returns:
        --------
        periodic_match : ndarray (n,)
            Indices of the coordinates which are identical by
            periodic boundary conditions.
        """
        if screen:
            original_index = np.arange(self.frac_coords.shape[0])[self.screen]
            coords = self.frac_coords[self.screen]
        else:
            original_index = np.arange(self.frac_coords.shape[0])
            coords = self.frac_coords
        periodic_match = original_index.copy()
        for i, j in enumerate(periodic_match):
            ind = original_index[i]
            if ind != j:
                continue

            new_match = matching_sites(self.frac_coords[j], coords)
            periodic_match[new_match] = ind

        return periodic_match

    def get_symmetric_sites(self, unique=True, screen=True):
        """Determine the symmetrically unique adsorption sites
        from a list of fractional coordinates.

        Parameters:
        -----------
        unique : bool
            Return only the unique symmetrically reduced sites.
        screen : bool
            Return only sites inside the unit cell.

        Returns:
        --------
        sites : dict of lists
            Dictionary of sites containing index of site
        """
        if screen is False:
            unique = False

        symmetry = utils.get_symmetry(self.slab, tol=self.tol)
        rotations = np.swapaxes(symmetry['rotations'], 1, 2)
        translations = symmetry['translations']
        affine = np.append(rotations, translations[:, None], axis=1)

        if screen:
            periodic = self.get_periodic_sites()
            points = self.frac_coords[periodic]
            true_index = np.where(self.screen)[0]
        else:
            points = self.frac_coords
            true_index = np.arange(self.frac_coords.shape[0])

        affine_points = np.insert(points, 3, 1, axis=1)
        operations = np.dot(affine_points, affine)
        symmetry_match = np.arange(points.shape[0])
        for i, j in enumerate(symmetry_match):
            if i != j:
                continue

            d = operations[i, :, None] - points
            d -= np.round(d)
            dind = np.where((np.abs(d) < self.tol).all(axis=2))[-1]
            symmetry_match[np.unique(dind)] = true_index[i]

        if unique:
            return np.unique(symmetry_match)
        else:
            return symmetry_match

    def get_adsorption_vectors(self):
        """Returns the vectors representing the furthest distance from
        the neighboring atoms.

        Returns:
        --------
        vectors : ndarray (n, 3)
            Adsorption vectors for surface sites.
        """
        top_coords = self.coordinates[self.connectivity == 1]
        vectors = np.empty((self.coordinates.shape[0], 3))

        for i, s in enumerate(self.coordinates):
            plane_points = np.array(
                list(self.r1_topology[i]) + list(self.r2_topology[i]),
                dtype=int)
            vectors[i] = utils.plane_normal(top_coords[plane_points])

        return vectors

    def get_adsorption_graph(self, unique=True):
        """Return the edges of adsorption sties defined as all regions
        with adjacent vertices.

        Parameters:
        -----------
        unique : bool
            Return only the unique edges.

        Returns:
        --------
        edges : ndarray (n, 2)
            All edges crossing ridge or vertices indexed by the expanded
            unit slab.
        """
        vt = Voronoi(self.coordinates[:, :2])
        regions = -np.ones((len(vt.regions) - 1, 12), dtype=int)
        for i, p in enumerate(vt.point_region):
            n = len(vt.regions[p])
            regions[i, :n] = vt.regions[p]

        site_id = self.get_symmetric_sites(unique=False, screen=False)
        site_id = site_id + self.connectivity / 10

        edges, uniques, sym_edge = [], [], []
        for i in np.where(self.screen)[0]:
            vert = np.unique(regions[i])[1:]

            for v in vert:
                neighbor_regions = np.where(regions == v)[0][:, None]
                new_edges = np.insert(neighbor_regions, 0, i, axis=1)

                for edge in new_edges:
                    edge = sorted(edge)
                    if edge in edges or edge[0] == edge[1]:
                        continue

                    sym = sorted(site_id[edge])
                    if sym in sym_edge:
                        uniques += [False]
                    else:
                        uniques += [True]
                        sym_edge += [sym]
                    edges += [edge]

        edges = np.array(edges)
        if unique:
            return edges[uniques]

        return edges

    def view(self):
        """Create a visualization of the sites."""
        top = self.connectivity == 1
        other = self.connectivity != 1
        dt = Delaunay(self.coordinates[:, :2][top])
        plt.triplot(dt.points[:, 0], dt.points[:, 1], dt.simplices.copy())
        plt.plot(dt.points[:, 0], dt.points[:, 1], 'o')
        plt.plot(self.coordinates[:, 0][other], self.coordinates[:, 1][other],
                 'o')
        plt.show()


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


def get_adsorption_sites(
        slab,
        symmetry_reduced=True,
        adsorption_vectors=False,
        tol=1e-5
):
    """Get the adsorption sites of a slab as defined by surface
    symmetries of the surface atoms.

    Parameters:
    -----------
    slab : atoms object
        The slab to find adsorption sites for. Must have surface
        atoms defined.
    symmetry_reduced : bool
        Return the symmetrically unique sites only.
    adsorption_vectors : bool
        Return the adsorption vectors.

    Returns:
    --------
    coordinates : ndarray (n, 3)
        Cartesian coordinates of activate sites.
    connectivity : ndarray (n,)
        Number of bonds formed for a given adsorbate.
    vectors : ndarray (n, 3)
        Vector associated with minimum surface interaction.
    """
    sites = AdsorptionSites(slab)

    if symmetry_reduced:
        s = sites.get_symmetric_sites()
    else:
        s = sites.get_periodic_sites()

    coordinates = sites.coordinates[s]
    connectivity = sites.connectivity[s]
    if adsorption_vectors:
        vectors = sites.get_adsorption_vectors()[s]

        return coordinates, connectivity, vectors

    return coordinates, connectivity
