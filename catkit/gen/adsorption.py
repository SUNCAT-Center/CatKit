from . import defaults
from . import utils
from . import symmetry
import catkit
import matplotlib.pyplot as plt
import itertools
import networkx as nx
import numpy as np
import scipy
radii = defaults.get('radii')


class AdsorptionSites():
    """Adsorption site object."""

    def __init__(self, slab, surface_atoms=None, tol=1e-5):
        """Create an extended unit cell of the surface sites for
        use in identifying other sites.

        Parameters
        ----------
        slab : Gatoms object
            The slab associated with the adsorption site network to be
            attached.
        tol : float
            Absolute tolerance for floating point errors.
        """
        index, coords, offsets = utils.expand_cell(slab, cutoff=5.0)
        if surface_atoms is None:
            surface_atoms = slab.get_surface_atoms()
        if surface_atoms is None:
            raise ValueError('Slab must contain surface atoms')

        extended_top = np.where(np.in1d(index, surface_atoms))[0]

        self.tol = tol
        self.coordinates = coords[extended_top].tolist()
        self.connectivity = np.ones(extended_top.shape[0]).tolist()
        self.r1_topology = [[i] for i in np.arange(len(extended_top))]
        self.index = index[extended_top]

        sites = self._get_higher_coordination_sites(coords[extended_top])
        self.r2_topology = sites['top'][2]

        # Put data into array format
        selection = ['bridge', 'hollow', '4fold']
        for i, k in enumerate(selection):
            coordinates, r1top, r2top = sites[k]

            if k in ['hollow', '4fold']:
                r2top = [[] for _ in coordinates]

            self.connectivity += (np.ones(len(coordinates)) * (i + 2)).tolist()
            self.coordinates += coordinates
            self.r1_topology += r1top
            self.r2_topology += r2top

        self.coordinates = np.array(self.coordinates)
        self.connectivity = np.array(self.connectivity, dtype=int)
        self.r1_topology = np.array(self.r1_topology, dtype=object)
        self.r2_topology = np.array(self.r2_topology, dtype=object)
        self.frac_coords = np.dot(self.coordinates, np.linalg.pinv(slab.cell))
        self.slab = slab

        screen = (self.frac_coords[:, 0] > 0 - self.tol) & \
                 (self.frac_coords[:, 0] < 1 - self.tol) & \
                 (self.frac_coords[:, 1] > 0 - self.tol) & \
                 (self.frac_coords[:, 1] < 1 - self.tol)

        self.screen = screen
        self._symmetric_sites = None

    def get_connectivity(self, unique=True):
        """Return the number of connections associated with each site."""
        if unique:
            sel = self.get_symmetric_sites()
        else:
            sel = self.get_periodic_sites()

        return self.connectivity[sel]

    def get_coordinates(self, unique=True):
        """Return the 3D coordinates associated with each site."""
        if unique:
            sel = self.get_symmetric_sites()
        else:
            sel = self.get_periodic_sites()

        return self.coordinates[sel]

    def get_topology(self, unique=True):
        """Return the indices of adjacent surface atoms."""
        topology = [self.index[top] for top in self.r1_topology]
        topology = np.array(topology)
        if unique:
            sel = self.get_symmetric_sites()
        else:
            sel = self.get_periodic_sites()

        return topology[sel]

    def _get_higher_coordination_sites(self,
                                       top_coordinates,
                                       allow_obtuse=True):
        """Find all bridge and hollow sites (3-fold and 4-fold) given an
        input slab based Delaunay triangulation of surface atoms of a
        super-cell.

        TODO: Determine if this can be made more efficient by
        removing the 'sites' dictionary.

        Parameters
        ----------
        top_coordinates : ndarray (n, 3)
            Cartesian coordinates for the top atoms of the unit cell.

        Returns
        -------
        sites : dict of 3 lists
            Dictionary sites containing positions, points, and neighbor lists.
        """
        sites = {
            'top': [top_coordinates, [], [[] for _ in top_coordinates]],
            'bridge': [[], [], []],
            'hollow': [[], [], []],
            '4fold': [[], [], []]}

        dt = scipy.spatial.Delaunay(sites['top'][0][:, :2])
        neighbors = dt.neighbors
        simplices = dt.simplices

        for i, corners in enumerate(simplices):
            cir = scipy.linalg.circulant(corners)
            edges = cir[:, 1:]

            # Inner angle of each triangle corner
            vec = sites['top'][0][edges.T] - sites['top'][0][corners]
            uvec = vec.T / np.linalg.norm(vec, axis=2).T
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
                    ouvec = ovec / np.linalg.norm(ovec)
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

            edges = itertools.product(s[:2], s[2:])
            for edge in edges:
                edge = sorted(edge)
                i = sites['bridge'][1].index(edge)
                n, m = sites['bridge'][1][i], sites['bridge'][2][i]
                nn = list(set(s) - set(n + m))

                if len(nn) == 0:
                    continue
                sites['bridge'][2][i] += [nn[0]]

        return sites

    def get_periodic_sites(self, screen=True):
        """Return an index of the coordinates which are unique by
        periodic boundary conditions.

        Parameters
        ----------
        screen : bool
            Return only sites inside the unit cell.

        Returns
        -------
        periodic_match : ndarray (n,)
            Indices of the coordinates which are identical by
            periodic boundary conditions.
        """
        periodic_match = np.arange(self.frac_coords.shape[0])

        if screen:
            return periodic_match[self.screen]

        coords = self.frac_coords
        periodic = periodic_match.copy()[self.screen]

        for p in periodic:
            matched = utils.matching_sites(self.frac_coords[p], coords)
            periodic_match[matched] = p

        return periodic_match

    def get_symmetric_sites(self, unique=True, screen=True):
        """Determine the symmetrically unique adsorption sites
        from a list of fractional coordinates.

        Parameters
        ----------
        unique : bool
            Return only the unique symmetrically reduced sites.
        screen : bool
            Return only sites inside the unit cell.

        Returns
        -------
        sites : dict of lists
            Dictionary of sites containing index of site
        """
        symmetry_match = self._symmetric_sites

        if symmetry_match is None:
            sym = symmetry.Symmetry(self.slab, tol=self.tol)

            rotations, translations = sym.get_symmetry_operations(affine=False)
            rotations = np.swapaxes(rotations, 1, 2)

            affine = np.append(rotations, translations[:, None], axis=1)
            points = self.frac_coords
            true_index = self.get_periodic_sites(False)

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

            self._symmetric_sites = symmetry_match

        if screen:
            periodic = self.get_periodic_sites()
            symmetry_match = symmetry_match[periodic]
        if unique:
            return np.unique(symmetry_match)

        else:
            return symmetry_match

    def get_adsorption_vectors(self, unique=True, screen=True):
        """Returns the vectors representing the furthest distance from
        the neighboring atoms.

        Returns
        -------
        vectors : ndarray (n, 3)
            Adsorption vectors for surface sites.
        """
        top_coords = self.coordinates[self.connectivity == 1]
        if unique:
            sel = self.get_symmetric_sites(screen=screen)
        else:
            sel = self.get_periodic_sites(screen=screen)
        coords = self.coordinates[sel]
        r1top = self.r1_topology[sel]
        r2top = self.r2_topology[sel]

        vectors = np.empty((coords.shape[0], 3))
        for i, s in enumerate(coords):
            plane_points = np.array(list(r1top[i]) + list(r2top[i]), dtype=int)
            vectors[i] = utils.plane_normal(top_coords[plane_points])

        return vectors

    def get_adsorption_edges(self, symmetric=True, periodic=True):
        """Return the edges of adsorption sties defined as all regions
        with adjacent vertices.

        Parameters
        ----------
        symmetric : bool
            Return only the symmetrically reduced edges.
        periodic : bool
            Return edges which are unique via periodicity.

        Returns
        -------
        edges : ndarray (n, 2)
            All edges crossing ridge or vertices indexed by the expanded
            unit slab.
        """
        vt = scipy.spatial.Voronoi(self.coordinates[:, :2],
                                   qhull_options='Qbb Qc Qz C{}'.format(1e-2))

        select, lens = [], []
        for i, p in enumerate(vt.point_region):
            select += [vt.regions[p]]
            lens += [len(vt.regions[p])]

        dmax = max(lens)
        regions = np.zeros((len(select), dmax), int)
        mask = np.arange(dmax) < np.array(lens)[:, None]
        regions[mask] = np.concatenate(select)

        site_id = self.get_symmetric_sites(unique=False, screen=False)
        site_id = site_id + self.connectivity / 10
        per = self.get_periodic_sites(screen=False)

        uper = self.get_periodic_sites()
        edges, symmetry, uniques = [], [], []
        for i, p in enumerate(uper):
            poi = vt.point_region[p]
            voi = vt.regions[poi]

            for v in voi:
                nr = np.where(regions == v)[0]

                for n in nr:
                    edge = sorted((p, n))

                    if n in uper[:i + 1] or edge in edges:
                        continue

                    if (np.in1d(per[edge], per[uper[:i]]).any()) and periodic:
                        continue

                    sym = sorted(site_id[edge])
                    if sym in symmetry:
                        uniques += [False]
                    else:
                        uniques += [True]
                        symmetry += [sym]

                    edges += [edge]

        edges = np.array(edges)
        if symmetric:
            edges = edges[uniques]

        return edges

    def ex_sites(self, index, select='inner', cutoff=0):
        """Get site indices inside or outside of a cutoff radii from a
        provided periodic site index. If two sites are provided, an
        option to return the mutually inclusive points is also available.
        """
        per = self.get_periodic_sites(False)
        sym = self.get_symmetric_sites()
        edges = self.get_adsorption_edges(symmetric=False, periodic=False)
        coords = self.coordinates[:, :2]

        if isinstance(index, int):
            index = [index]

        if not cutoff:
            for i in per[index]:
                sd = np.where(edges == i)[0]
                select_coords = coords[edges[sd]]
                d = np.linalg.norm(np.diff(select_coords, axis=1), axis=2)
                cutoff = max(d.max(), cutoff)
        cutoff += self.tol

        diff = coords[:, None] - coords[sym]
        norm = np.linalg.norm(diff, axis=2)
        neighbors = np.array(np.where(norm < cutoff))

        neighbors = []
        for i in index:
            diff = coords[:, None] - coords[per[i]]
            norm = np.linalg.norm(diff, axis=2)
            if select == 'mutual' and len(index) == 2:
                neighbors += [np.where(norm < cutoff)[0].tolist()]
            else:
                neighbors += np.where(norm < cutoff)[0].tolist()

        if select == 'inner':
            return per[neighbors]
        elif select == 'outer':
            return np.setdiff1d(per, per[neighbors])
        elif select == 'mutual':
            return np.intersect1d(per[neighbors[0]], per[neighbors[1]])

    def plot(self, savefile=None):
        """Create a visualization of the sites."""
        top = self.connectivity == 1
        other = self.connectivity != 1
        dt = scipy.spatial.Delaunay(self.coordinates[:, :2][top])

        fig = plt.figure(figsize=(6, 4), frameon=False)
        ax = fig.add_axes([0, 0, 1, 1])

        ax.triplot(dt.points[:, 0], dt.points[:, 1], dt.simplices.copy())
        ax.plot(dt.points[:, 0], dt.points[:, 1], 'o')
        ax.plot(self.coordinates[:, 0][other], self.coordinates[:, 1][other],
                'o')
        ax.axis('off')

        if savefile:
            plt.savefig(savefile, transparent=True)
        else:
            plt.show()


class Builder(AdsorptionSites):
    """Initial module for creation of 3D slab structures with
    attached adsorbates.
    """

    def __repr__(self):
        formula = self.slab.get_chemical_formula()
        string = 'Adsorption builder for {} slab.\n'.format(formula)
        sym = len(self.get_symmetric_sites())
        string += 'unique adsorption sites: {}\n'.format(sym)
        con = self.get_connectivity()
        string += 'site connectivity: {}\n'.format(con)
        edges = self.get_adsorption_edges()
        string += 'unique adsorption edges: {}'.format(len(edges))

        return string

    def add_adsorbates(self, adsorbates, indices):
        """ """
        slab = self.slab.copy()
        for i, ads in enumerate(adsorbates):
            bond = np.where(ads.get_tags() == -1)[0][0]

            slab = self._single_adsorption(
                ads,
                slab=slab,
                bond=bond,
                site_index=indices[i],
                symmetric=False)

        return slab

    def add_adsorbate(
            self,
            adsorbate,
            bonds=None,
            index=0,
            auto_construct=True,
            **kwargs):
        """Add and adsorbate to a slab. If the auto_constructor flag is False,
        the atoms object provided will be attached at the active site.

        Parameters
        ----------
        adsorbate : gratoms object
            Molecule to connect to the surface.
        bonds : int or list of 2 int
            Index of adsorbate atoms to be bonded.
        index : int
            Index of the site or edge to use as the adsorption position. A
            value of -1 will return all possible structures.
        auto_construct : bool
            Whether to automatically estimate the position of atoms in larger
            molecules or use the provided structure.

        Returns
        -------
        slabs : gratoms object
            Slab(s) with adsorbate attached.
        """
        if bonds is None:
            # Molecules with tag -1 are designated to bond
            bonds = np.where(adsorbate.get_tags() == -1)[0]

        if len(bonds) == 0:
            raise ValueError('Specify the index of atom to bond.')

        elif len(bonds) == 1:
            if index == -1:
                slab = []
                for i, _ in enumerate(self.get_symmetric_sites()):
                    slab += [self._single_adsorption(
                        adsorbate,
                        bond=bonds[0],
                        site_index=i,
                        auto_construct=auto_construct,
                        **kwargs)]
            elif isinstance(index, (list, np.ndarray)):
                slab = []
                for i in index:
                    slab += [self._single_adsorption(
                        adsorbate,
                        bond=bonds[0],
                        site_index=i,
                        auto_construct=auto_construct,
                        **kwargs)]
            else:
                slab = self._single_adsorption(
                    adsorbate,
                    bond=bonds[0],
                    site_index=index,
                    auto_construct=auto_construct,
                    **kwargs)

        elif len(bonds) == 2:
            if index == -1:
                slab = []
                edges = self.get_adsorption_edges()
                for i, _ in enumerate(edges):
                    slab += [self._double_adsorption(
                        adsorbate,
                        bonds=bonds,
                        edge_index=i,
                        **kwargs)]
            else:
                slab = self._double_adsorption(
                    adsorbate,
                    bonds=bonds,
                    edge_index=index,
                    **kwargs)

        else:
            raise ValueError('Only mono- and bidentate adsorption supported.')

        return slab

    def _single_adsorption(
            self,
            adsorbate,
            bond,
            slab=None,
            site_index=0,
            auto_construct=True,
            symmetric=True):
        """Bond and adsorbate by a single atom."""
        if slab is None:
            slab = self.slab.copy()
        atoms = adsorbate.copy()
        atoms.set_cell(slab.cell)

        if symmetric:
            ind = self.get_symmetric_sites()[site_index]
            vector = self.get_adsorption_vectors()[site_index]
        else:
            ind = self.get_periodic_sites()[site_index]
            vector = self.get_adsorption_vectors(unique=False)[site_index]

        # Improved position estimate for site.
        u = self.r1_topology[ind]
        r = radii[slab[self.index[u]].numbers]
        top_sites = self.coordinates[self.connectivity == 1]

        numbers = atoms.numbers[bond]
        R = radii[numbers]
        base_position = utils.trilaterate(top_sites[u], r + R, vector)

        branches = nx.bfs_successors(atoms.graph, bond)
        atoms.translate(-atoms.positions[bond])

        if auto_construct:
            atoms = catkit.gen.molecules.get_3D_positions(atoms, bond)

            # Align with the adsorption vector
            atoms.rotate([0, 0, 1], vector)

        atoms.translate(base_position)
        n = len(slab)
        slab += atoms

        # Add graph connections
        for metal_index in self.index[u]:
            slab.graph.add_edge(metal_index, bond + n)

        return slab

    def _double_adsorption(self, adsorbate,
                           bonds=None, edge_index=0):
        """Bond and adsorbate by two adjacent atoms."""
        slab = self.slab.copy()
        atoms = adsorbate.copy()
        atoms.set_cell(slab.cell)

        numbers = atoms.numbers[bonds]
        R = radii[numbers] * 0.95
        edges = self.get_adsorption_edges()
        coords = self.coordinates[edges[edge_index]]

        U = self.r1_topology[edges[edge_index]]
        for i, u in enumerate(U):
            r = radii[slab[self.index[u]].numbers] * 0.95
            top_sites = self.coordinates[self.connectivity == 1]
            coords[i] = utils.trilaterate(top_sites[u], R[i] + r)

        vec = coords[1] - coords[0]
        n = np.linalg.norm(vec)
        uvec0 = vec / n
        d = np.sum(radii[numbers]) * 0.95
        dn = (d - n) / 2

        base_position0 = coords[0] - uvec0 * dn
        base_position1 = coords[1] + uvec0 * dn

        # Position the base atoms
        atoms[bonds[0]].position = base_position0
        atoms[bonds[1]].position = base_position1

        # Temporarily break adsorbate bond
        atoms.graph.remove_edge(*bonds)

        vectors = self.get_adsorption_vectors(screen=False, unique=False)
        uvec1 = vectors[edges[edge_index]]
        uvec2 = np.cross(uvec1, uvec0)
        uvec2 /= -np.linalg.norm(uvec2, axis=1)[:, None]
        uvec1 = np.cross(uvec2, uvec0)

        branches0 = list(nx.bfs_successors(atoms.graph, bonds[0]))
        if len(branches0[0][1]) != 0:
            uvec = [-uvec0, uvec1[0], uvec2[0]]
            self._branch_bidentate(atoms, uvec, branches0[0])
            for branch in branches0[1:]:
                catkit.gen.molecules._branch_molecule(
                    atoms, branch, adsorption=True)

        branches1 = list(nx.bfs_successors(atoms.graph, bonds[1]))
        if len(branches1[0][1]) != 0:
            uvec = [uvec0, uvec1[0], uvec2[0]]
            self._branch_bidentate(atoms, uvec, branches1[0])
            for branch in branches1[1:]:
                catkit.gen.molecules._branch_molecule(
                    atoms, branch, adsorption=True)

        n = len(slab)
        slab += atoms
        # Add graph connections
        slab.graph.add_edge(*np.array(bonds) + n)
        for i, u in enumerate(U):
            for metal_index in self.index[u]:
                slab.graph.add_edge(metal_index, bonds[i] + n)

        return slab

    def _branch_bidentate(self, atoms, uvec, branch):
        """Return extended positions for additional adsorbates
        based on provided unit vectors.
        """
        r, nodes = branch
        num = atoms.numbers[[r] + nodes]
        d = radii[num[1:]] + radii[num[0]]
        c = atoms[r].position

        # Single additional atom
        if len(nodes) == 1:
            coord0 = c + \
                d[0] * uvec[0] * np.cos(1 / 3. * np.pi) + \
                d[0] * uvec[1] * np.sin(1 / 3. * np.pi)
            atoms[nodes[0]].position = coord0

        # Two branch system
        elif len(nodes) == 2:
            coord0 = c + \
                d[0] * uvec[1] * np.cos(1 / 3. * np.pi) + \
                0.866 * d[0] * uvec[0] * np.cos(1 / 3. * np.pi) + \
                0.866 * d[0] * uvec[2] * np.sin(1 / 3. * np.pi)
            atoms[nodes[0]].position = coord0

            coord1 = c + \
                d[1] * uvec[1] * np.cos(1 / 3. * np.pi) + \
                0.866 * d[1] * uvec[0] * np.cos(1 / 3. * np.pi) + \
                0.866 * d[1] * -uvec[2] * np.sin(1 / 3. * np.pi)
            atoms[nodes[1]].position = coord1

        else:
            raise ValueError('Too many bonded atoms to position correctly.')


def get_adsorption_sites(slab,
                         surface_atoms=None,
                         symmetry_reduced=True,
                         tol=1e-5):
    """Get the adsorption sites of a slab as defined by surface
    symmetries of the surface atoms.

    Parameters
    ----------
    slab : Atoms object
        The slab to find adsorption sites for. Must have surface
        atoms defined.
    surface_atoms : array_like (N,)
        List of the surface atom indices.
    symmetry_reduced : bool
        Return the symmetrically unique sites only.
    adsorption_vectors : bool
        Return the adsorption vectors.

    Returns
    -------
    coordinates : ndarray (N, 3)
        Cartesian coordinates of activate sites.
    connectivity : ndarray (N,)
        Number of bonds formed for a given adsorbate.
    symmetry_index : ndarray (N,)
        Arbitrary indices of symmetric sites.
    """
    if surface_atoms is not None:
        slab.set_surface_atoms(surface_atoms)

    sites = AdsorptionSites(slab)

    if symmetry_reduced:
        s = sites.get_symmetric_sites()
    else:
        s = sites.get_periodic_sites()

    coordinates = sites.coordinates[s]
    connectivity = sites.connectivity[s]

    if symmetry_reduced:
        return coordinates, connectivity

    symmetries = sites.get_symmetric_sites(unique=False)
    unique, counts = np.unique(symmetries, return_counts=True)
    symmetry_index = np.arange(len(unique)).repeat(counts)

    return coordinates, connectivity, symmetry_index


def _get_adsorption_sites(surface_positions, tol=1e-5):
    """Return the positions and topology of adsorption sites based on the
    provided 3D positions of surface atoms for a structure.

    This function is intended for internal use only.

    Parameters
    ----------
    surface_positions : ndarray (N, 3)
        Cartesian coordinates for the surface atoms of a structure.
    tol : float
        Float point precision tolerance.

    Returns
    -------
    positions : ndarray (M, 3)
        Cartesian coordinates for the identifed adsorption sites. This is
        a superset which contrains the surface atoms positions.
    r1topology : list of lists (M, X)
        First nearest-neighbors of each of the adsorption sites associated
        with the positions.
    r2topology : list of lists (M, Y)
        Second nearest-neighbors of each of the adsorption sites associated
        with the positions.
    """
    positions = [surface_positions, [], [], []]
    r1topology = [[[i] for i, _ in enumerate(surface_positions)], [], [], []]
    r2topology = [[[] for _ in surface_positions], [], [], []]

    dt = scipy.spatial.Delaunay(positions[0][:, :2])
    neighbors = dt.neighbors
    simplices = dt.simplices

    for i, corners in enumerate(simplices):
        cir = scipy.linalg.circulant(corners)
        edges = cir[:, 1:]

        # Inner angle of each triangle corner
        vec = positions[0][edges.T] - positions[0][corners]
        uvec = vec.T / np.linalg.norm(vec, axis=2).T
        angles = np.sum(uvec.T[0] * uvec.T[1], axis=1)

        # Angle types
        right = np.isclose(angles, 0)
        obtuse = (angles < -tol)

        rh_corner = corners[right]
        edge_neighbors = neighbors[i]

        bridge = np.sum(positions[0][edges], axis=1) / 2

        # Looping through corners allows for elimination of
        # redundant points, identification of 4-fold hollows,
        # and collection of bridge neighbors.
        for j, c in enumerate(corners):
            edge = sorted(edges[j])

            if edge in r1topology[1]:
                continue

            # Get the bridge neighbors (for adsorption vector)
            neighbor_simplex = simplices[edge_neighbors[j]]
            oc = list(set(neighbor_simplex) - set(edge))[0]

            # Right angles potentially indicate 4-fold hollow
            potential_hollow = edge + sorted([c, oc])
            if c in rh_corner:

                if potential_hollow in r1topology[3]:
                    continue

                # Assumption: If not 4-fold, this suggests
                # no hollow OR bridge site is present.
                ovec = positions[0][edge] - positions[0][oc]
                ouvec = ovec / np.linalg.norm(ovec)
                oangle = np.dot(*ouvec)
                oright = np.isclose(oangle, 0)
                if oright:
                    positions[3] += [bridge[j]]
                    r1topology[3] += [potential_hollow]
                    r2topology[3] += []
                    r2topology[0][c] += [oc]
            else:
                positions[1] += [bridge[j]]
                r1topology[1] += [edge]
                r2topology[1] += [[c, oc]]

            r2topology[0][edge[0]] += [edge[1]]
            r2topology[0][edge[1]] += [edge[0]]

        if not right.any() and not obtuse.any():
            hollow = np.average(positions[0][corners], axis=0)
            positions[2] += [hollow]
            r1topology[2] += [corners.tolist()]
            r2topology[2] += []

    # For collecting missed bridge neighbors
    for s in r1topology[3]:

        edges = itertools.product(s[:2], s[2:])
        for edge in edges:
            edge = sorted(edge)
            i = r1topology[1].index(edge)
            n, m = r1topology[1][i], r2topology[1][i]
            nn = list(set(s) - set(n + m))

            if len(nn) == 0:
                continue
            r2topology[1][i] += [nn[0]]

    positions = np.array([j for i in positions for j in i])
    r1topology = np.array([j for i in r1topology for j in i])
    r2topology = np.array([j for i in r2topology for j in i])

    return positions, r1topology, r2topology


def symmetry_equivalent_points(
        fractional_coordinates, atoms, tol=1e-5):
    """Return the symmetrically equivalent points from a list of
    provided fractional coordinates.

    Parameters
    ----------
    fractional_coordinates : ndarray (N ,3)
        Fractional coordinates to find symmetrical equivalence
        between.
    atoms : Atoms object
        Atoms object to use the unit cell, positions, and pbc of.
    tol : float
        Float point precision tolerance.

    Returns
    -------
    symmetry_match : ndarray (N,)
        Indices of fractional coordinates which are unique or
        matching.
    """
    sym = symmetry.Symmetry(atoms, tol=1e-2)
    affine_matrices = sym.get_symmetry_operations()
    affine_matrices = np.swapaxes(affine_matrices, 1, 2)

    pbc = ~np.array(atoms.pbc.tolist() + [True])
    nt = np.isclose(affine_matrices[:, pbc, 3], 0).any(axis=1)
    nt &= np.isclose(affine_matrices[:, pbc, 2], 1).any(axis=1)
    affine_matrices = affine_matrices[nt]

    affine_points = np.insert(fractional_coordinates, 3, 1, axis=1)
    operations = np.dot(affine_points, affine_matrices)[:, :, :3]

    periodic_match = np.arange(fractional_coordinates.shape[0])
    symmetry_match = periodic_match.copy()
    for i, j in enumerate(symmetry_match):
        if i != j:
            continue

        d = operations[i, :, None] - fractional_coordinates
        d -= np.round(d)
        dind = np.where((np.abs(d) < tol).all(axis=2))[-1]
        symmetry_match[np.unique(dind)] = periodic_match[i]

    return symmetry_match
