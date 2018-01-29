from . import utils
import numpy as np
from numpy.linalg import norm, pinv, solve
from ase.neighborlist import NeighborList
import networkx as nx
import networkx.algorithms.isomorphism as iso
from ase.build import rotate
from ase.constraints import FixAtoms
from ase import Atoms
from scipy.spatial import Delaunay
from scipy.linalg import circulant
from itertools import product
try:
    from math import gcd
except ImportError:
    from fractions import gcd


class SlabGenerator(object):
    """ Class for generation of slab unit cells from bulk unit cells.
    """

    def __init__(
            self,
            bulk,
            miller_index=[1, 1, 1],
            layers=4,
            fixed=2,
            vacuum=0,
            tol=1e-8
    ):
        """ Generate a slab from an ASE bulk atoms-object.

        Parameters:
          bulk: ASE atoms-object
            Bulk structure to produce the slab from.
          miller_index: list (3,)
            Miller index to construct surface from.
          layers: int
            Number of layers to include in the slab.
          fixed: int
            Number of layers to fix in the slab.
          vacuum: float
            Angstroms of vacuum to add to the slab.
          tol: float
            Tolerance for floating point rounding errors.
        """

        self.bulk = bulk
        self.miller_index = np.array(miller_index)
        self.layers = layers
        self.fixed = fixed
        self.vacuum = vacuum
        self.tol = tol

        self.unique_terminations = None
        self.surface_atoms = None

        self._basis = self.build_basis()

    def build_basis(self):
        """ Get the basis unit cell from bulk unit cell. This
        basis is effectively the same as the bulk, but rotated such
        that the z-axis is aligned with the surface termination.

        The basis is stored separately from the slab generated. This
        is temporary until a slab class is created.

        Returns:
          basis: ASE atoms-object
            The basis slab corresponding to the provided bulk.
        """

        h, k, l = self.miller_index
        h0, k0, l0 = (self.miller_index == 0)
        if h0 and k0 or h0 and l0 or k0 and l0:
            if not h0:
                c1, c2, c3 = [(0, 1, 0), (0, 0, 1), (1, 0, 0)]
            if not k0:
                c1, c2, c3 = [(0, 0, 1), (1, 0, 0), (0, 1, 0)]
            if not l0:
                c1, c2, c3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        else:
            p, q = ext_gcd(k, l)
            a1, a2, a3 = self.bulk.cell

            # constants describing the dot product of basis c1 and c2:
            # dot(c1,c2) = k1+i*k2, i in Z
            k1 = np.dot(p * (k * a1 - h * a2) + q * (l * a1 - h * a3),
                        l * a2 - k * a3)
            k2 = np.dot(l * (k * a1 - h * a2) - k * (l * a1 - h * a3),
                        l * a2 - k * a3)

            if abs(k2) > self.tol:
                # i corresponding to the optimal basis
                i = -int(round(k1 / k2))
                p, q = p + i * l, q - i * k

            a, b = ext_gcd(p * k + q * l, h)

            c1 = (p * k + q * l, -p * h, -q * h)
            c2 = np.array((0, l, -k)) // abs(gcd(l, k))
            c3 = (b, a * p, a * q)

        slab = self.bulk.copy()
        basis = np.array([c1, c2, c3])

        scaled = solve(basis.T, slab.get_scaled_positions().T).T
        scaled -= np.floor(scaled + self.tol)
        slab.set_scaled_positions(scaled)
        slab.set_cell(np.dot(basis, slab.cell), scale_atoms=True)

        a1, a2, a3 = slab.cell
        a3 = np.cross(a1, a2) / norm(np.cross(a1, a2))
        rotate(slab, a3, (0, 0, 1), a1, (1, 0, 0))

        return slab

    def get_unique_terminations(self):
        """ Return smallest unit cell corresponding to given surface and
        unique surface terminations based on symmetry and nearest neighbors.

        Returns:
          unique_terminations: list
            Unique terminations of a surface.
        """

        # Find all different planes as simply different z-coordinates
        z_planes = utils.get_unique_coordinates(self._basis, tol=self.tol)

        # now get the symmetries of lattice
        symmetry = utils.get_symmetry(self._basis, tol=self.tol)
        rotations = symmetry['rotations']
        translations = symmetry['translations']

        # Find all symmetries which are rotations about the z-axis
        z_symmetry = []
        for i, rotation in enumerate(rotations):
            if (abs(rotation[2][0]) < self.tol and
                abs(rotation[2][1]) < self.tol and
                abs(rotation[0][2]) < self.tol and
                abs(rotation[1][2]) < self.tol and
                abs(rotation[2][2] - 1.0) < self.tol):

                if not np.isclose(
                        translations[i][2],
                        z_symmetry,
                        rtol=self.tol).any():
                    z_symmetry += [translations[i][2]]

        # Find all unique z-shifts
        unique_shift = [z_planes[0]]
        for i in range(1, len(z_planes)):
            symmetry_found = False
            for j in range(0, i):
                z_diff = z_planes[i] - z_planes[j]
                for z_sym in z_symmetry:
                    if np.allclose(z_sym, z_diff, rtol=self.tol):
                        symmetry_found = True
                        break
                else:
                    continue
                break

            if not symmetry_found:
                unique_shift += [z_planes[i]]

        if len(unique_shift) == 1:
            return unique_shift

        # Now search symmetrically unique planes
        # For nearest-neighbor uniqueness
        unique_terminations, graphs = [], []
        for i, z_shift in enumerate(unique_shift):
            tmp_slab = self._basis.copy()
            tmp_slab.translate([0, 0, -z_shift])
            tmp_slab.wrap(pbc=[1, 1, 1])

            zpos = tmp_slab.get_scaled_positions()[:, 2]
            index = np.arange(len(tmp_slab))
            del tmp_slab[index[zpos < 0.5]]

            nl = NeighborList(
                [2] * len(tmp_slab),
                skin=0.0,
                bothways=True,
                self_interaction=False)
            nl.build(tmp_slab)

            G = nx.MultiGraph()
            symbols = tmp_slab.get_chemical_symbols()
            for node, neighbors in enumerate(nl.neighbors):
                G.add_node(node, symbols=symbols[node])
                d = tmp_slab.get_distances(node, neighbors, mic=True)
                edges = [[node, _, {'distance': d[i]}] for i, _ in
                         enumerate(nl.get_neighbors(node)[0])]
                G.add_edges_from(edges)

            isomorph = False
            for G0 in graphs:
                nm = iso.categorical_node_match('symbols', 'X')
                em = iso.numerical_multiedge_match('distance', 1)
                if nx.is_isomorphic(
                        G, G0,
                        edge_match=em,
                        node_match=nm):
                    isomorph = True
                    break

            if not isomorph:
                graphs += [G]
                unique_terminations += [z_shift]

        self.unique_terminations = unique_terminations

        return unique_terminations

    def get_slab(self, iterm=None, primitive=False):
        """ Generate a slab object with a certain number of layers.

        Parameters:
          primitive: bool
            Whether to reduce the unit cell to its primitive form.
          iterm: int
            A termination index in reference to the list of possible
            terminations.

        Returns: ASE atoms-object
          The modified basis slab produced based on the layer specifications
          given.
        """

        slab = self._basis.copy()

        if iterm:
            if self.unique_terminations is None:
                terminations = self.get_unique_terminations()
            else:
                terminations = self.unique_terminations
            zshift = terminations[iterm]

            slab.translate([0, 0, -zshift])
            slab.wrap(pbc=True)

        # Get the minimum number of layers needed
        zlayers = utils.get_unique_coordinates(
            slab,
            direct=False,
            tol=self.tol
        )
        z_repetitions = np.ceil(self.layers / len(zlayers))
        slab *= (1, 1, int(z_repetitions))

        # Orthogonalize the z-coordinate
        # Warning: bulk symmetry is lost at this point
        a1, a2, a3 = slab.cell
        a3 = (np.cross(a1, a2) * np.dot(a3, np.cross(a1, a2)) /
              norm(np.cross(a1, a2)) ** 2)
        slab.cell[2] = a3

        if primitive:
            if self.vacuum:

                slab.center(vacuum=self.vacuum, axis=2)
            else:
                raise(
                    NotImplementedError,
                    'Primitive slab generation requires vacuum'
                )

            slab = utils.get_primitive_cell(slab)

            # For hcp(1, 1, 0), primitive alters z-axis
            d = norm(slab.cell, axis=0)
            maxd = np.argwhere(d == d.max())[0][0]
            if maxd != 2:
                slab.rotate(slab.cell[maxd], 'z', rotate_cell=True)
                slab.cell[[maxd, 2]] = slab.cell[[2, maxd]]
                slab.cell[maxd] = -slab.cell[maxd]
                slab.wrap(pbc=True)

            slab.rotate(slab.cell[0], 'x', rotate_cell=True)

            # spglib occasionally returns a bimodal slab
            zpos = slab.get_scaled_positions().T[2]
            if zpos.max() > 0.9 or zpos.min() < 0.1:
                translate = slab.positions.T[2][zpos > 0.5].min()
                slab.positions -= [0, 0, translate + self.tol]
                slab.wrap(pbc=True)
                slab.center(vacuum=self.vacuum, axis=2)

        # Get the direct z-coordinate of the requested layer
        zlayers = utils.get_unique_coordinates(
            slab,
            direct=False,
            tag=True,
            tol=self.tol
        )
        ncut = sorted(zlayers)[::-1][:self.layers][-1]

        zpos = slab.positions[:, 2]
        index = np.arange(len(slab))
        del slab[index[zpos - ncut < -self.tol]]

        slab.cell[2][2] -= ncut
        slab.translate([0, 0, -ncut])

        del slab.constraints

        tags = slab.get_tags()
        fix = tags.max() - self.fixed

        constraints = FixAtoms(indices=[a.index for a in slab if a.tag > fix])
        slab.set_constraint(constraints)

        if self.vacuum:
            slab.center(vacuum=self.vacuum, axis=2)

        slab.wrap()
        slab.pbc = [1, 1, 0]

        return slab

    def get_surface_atoms(self, slab=None):
        """ Find the under-coordinated atoms at the upper and lower
        fraction of a given unit cell based on the bulk structure it
        originated from.

        Assumes the xy-plane is perpendicular to the miller index.

        Parameters:
          slab: ASE atoms-object
            The slab to find top layer atoms from.

        Returns: ndarray (n,), ndarray (m,)
          Array of atom indices corresponding to the top and bottom
          layers of the slab.
        """

        if slab is None:
            slab = self.get_slab(primitive=True)

        ind, N, mcut = utils.get_voronoi_neighbors(self.bulk)
        ind0, N0 = utils.get_cutoff_neighbors(slab, cutoff=mcut)

        ind = np.repeat(ind, np.ceil(len(ind0) / len(ind)))

        surf_atoms = np.nonzero(ind0 - ind[:len(ind0)])[0]

        hwp = slab.positions[surf_atoms] - slab.get_center_of_mass()
        top = surf_atoms[hwp.T[2] > 0]
        bottom = surf_atoms[hwp.T[2] < 0]

        self.surface_atoms = top

        return top, bottom

    def get_adsorption_sites(
            self,
            slab=None,
            surface_sites=None,
            symmetry_reduced=True,
            vectors=False,
            return_proxy=False,
    ):
        """ Helper function for getting the adsorption sites of a
        slab as defined by surface atom symmetries.

        Parameters:
          slab: ASE atoms-object
            The slab to find adsorption sites for.
          surface_sites: ndarray (n,)
            surface sites of the provided slab.
          symmetry_reduced: int
            Whether to return the symmetrically unique sites only.
          vectors: bool
            Whether to compute the adsorption vectors.
          return_proxy: bool
            Whether to return the proxy_slab for adsorption work.

        Returns: dict of 3 lists
          Dictionary of top, bridge, hollow, and 4-fold sites containing
          positions, points, and neighbor lists. If adsorption vectors
          are requested, the third list is replaced.
        """

        if slab is None:
            slab = self.get_slab(primitive=True)

        if surface_sites is None:
            if self.surface_atoms is None:
                surface_sites, _ = self.get_surface_atoms(slab)

        sites, vslab = find_adsorption_sites(
            slab=slab,
            surface_sites=surface_sites,
            trim=0.5,
            tol=self.tol,
        )

        sites = get_reduced_sites(
            sites=sites,
            slab=slab,
            tol=self.tol,
        )

        if symmetry_reduced:
            sites = get_symmetric_sites(
                sites=sites,
                slab=slab,
                tol=self.tol,
            )

        if vectors:
            sites = self._get_adsorption_vectors(vslab, sites)

        if return_proxy:
            return sites, vslab

        return sites

    def _get_adsorption_vectors(
            self,
            vslab,
            sites,
    ):
        """ Returns the vectors representing the furthest distance from
        the neighboring atoms.

        (TODO: This input is complex and confusing. Would be nice to simplify.)

        Parameters:
          vslab: ASE atoms-object
            The virtual surface produced from find_adsorption_sites.
          sites: dict of 3 lists
            Dictionary of top, bridge, hollow, and 4-fold sites containing
            positions, points, and neighbor lists.

        Returns: dict of 3 lists
          Dictionary of top, bridge, hollow, and 4-fold sites containing
          positions, points, and adsorption vector lists.
        """

        pos = vslab.positions
        for k, v in sites.items():
            coordinates, points, neighbors = v

            vectors = []
            for i, s in enumerate(coordinates):

                if len(neighbors):
                    direct = pos[np.append(points[i], neighbors[i])]
                else:
                    direct = pos[points[i]]

                vectors += [utils.plane_normal(direct)]

            sites[k][2] = np.array(vectors)

        return sites


def find_adsorption_sites(
        slab,
        surface_sites,
        trim=0.5,
        tol=1e-5,
):
    """ Find all bridge and hollow sites (3-fold and 4-fold) given an
    input slab based Delaunay triangulation of surface atoms of a
    super-cell.

    Parameters:
      slab: ASE atoms-object
        The slab to find adsorption sites for.
      surface_sites: ndarray (n,)
        surface sites of the provided slab.
      trim: float
        Percentage of fractional coordinates to remove.
      tol: float
        Absolute tolerance for floating point errors.

    Returns: dict of 3 lists
      Dictionary of top, bridge, hollow, and 4-fold sites containing
      positions, points, and neighbors of the super-cell.
    """

    # Top sites projected into expanded unit cell
    index, coords, offsets = utils.expand_cell(slab, r=8)

    # Some magic to get the top sites indices
    n = len(slab)
    rtop = surface_sites[:, None].repeat(len(offsets), axis=1)
    rtop = (rtop + np.arange(0, len(coords), n)).T.flatten()

    top_coords = coords[rtop]
    numbers = slab.get_atomic_numbers().tolist() * len(offsets)

    vslab = Atoms(
        positions=top_coords,
        symbols=np.array(numbers)[rtop],
        cell=slab.cell,
        tags=rtop,
    )

    # Dict of 2 lists (positions, site_id)
    sites = {
        'top': [
            top_coords,
            np.arange(rtop.shape[0]),
            [[] for _ in range(rtop.shape[0])]
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
        obtuse = (angles < -tol)

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
            # print(potential_hollow)
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

    xlim = offsets.T[0].max()
    ylim = offsets.T[1].max()

    # Convert lists to arrays
    for k, v in sites.items():
        positions, points, _ = v

        if len(positions) == 0:
            continue

        frac_coords = np.dot(positions, pinv(slab.cell))

        screen = (
            frac_coords[:, 0] < xlim + 1 - (trim * xlim)) & \
            (frac_coords[:, 0] > -xlim + (trim * xlim)) & \
            (frac_coords[:, 1] < ylim + 1 - (trim * ylim)) & \
            (frac_coords[:, 1] > -ylim + (trim * ylim))

        sites[k][0] = np.array(positions)[screen]
        sites[k][1] = np.array(points)[screen]

        if k in ['top', 'bridge']:
            sites[k][2] = np.array(sites[k][2])[screen]

    return sites, vslab


def get_reduced_sites(
        sites,
        slab,
        tol=1e-5,
):
    """ Reduce overlapping points via fractional coordinates. Intended
    for use after finding super-cell sites.

    Parameters:
      sites: dict
        Dictionary of cartesian coordinates to provide reduced sites from.
        Must have the form {'site': [[positions], [points], [neighbors]]}.
      slab: ASE atoms-object
        slab to determine symmetry operations from.
      tol: float
        Absolute tolerance for floating point errors.

    Returns: dict of 3 lists
      Dictionary sites containing positions, points, and neighbor lists.
    """

    for k, v in sites.items():
        positions, points, nnneighbors = v

        if len(positions) == 0:
            continue

        frac_coords = np.dot(positions, pinv(slab.cell))

        non_unique, unique = [], []
        for i, xyz in enumerate(frac_coords):
            if i in non_unique:
                continue

            non_unique += matching_sites(xyz, frac_coords).tolist()
            unique += [i]

        unique = np.array(unique)
        unique_positions = frac_coords[unique]

        shift = [tol, tol, 0]
        unique_positions += shift
        unique_positions %= 1
        unique_positions -= shift

        srt = np.lexsort((
            unique_positions[:, 0],
            unique_positions[:, 1]
        ))

        sites[k][0] = np.dot(unique_positions, slab.cell)[srt]
        sites[k][1] = points[unique][srt]

        if isinstance(nnneighbors, np.ndarray):
            sites[k][2] = nnneighbors[unique][srt]

    return sites


def get_symmetric_sites(
        sites,
        slab,
        tol=1e-5
):
    """ Determine the symmetrically unique adsorption sites
    from a dictionary of possible sites for a given slab.

    Parameters:
      sites: dict
        Dictionary of Cartesian coordinates to provide reduced sites from.
        Must have the form {'site': [[positions], [points], [neighbors]]}.
      slab: ASE atoms-object
        slab to determine symmetry operations from.
      tol: float
        Absolute tolerance for floating point errors.

    Returns: dict of 3 lists
      Dictionary sites containing positions, points, and neighbor lists.
    """

    symmetry = utils.get_symmetry(slab, tol=tol)
    rotations = symmetry['rotations']
    translations = symmetry['translations']

    for k, v in sites.items():
        positions, points, nnneighbors = v

        if len(positions) == 0:
            continue

        frac_coords = np.dot(positions, pinv(slab.cell))

        # Convert to fractional
        unique_positions, unique = [], []
        for i, xyz in enumerate(frac_coords):

            symmetry_match = False
            for j, rotation in enumerate(rotations):
                translation = translations[j]

                affine_matrix = np.eye(4)
                affine_matrix[0:3][:, 0:3] = rotation
                affine_matrix[0:3][:, 3] = translation

                affine_point = np.array([xyz[0], xyz[1], xyz[2], 1])
                operation = np.dot(affine_matrix, affine_point)[0:3]

                if len(matching_sites(operation, unique_positions)) > 0:
                    symmetry_match = True
                    break

            if not symmetry_match:
                unique_positions += [xyz]
                unique += [i]

        unique = np.array(unique)
        sites[k][0] = np.dot(unique_positions, slab.cell)
        sites[k][1] = np.array(points[unique])

        if isinstance(nnneighbors, np.ndarray):
            sites[k][2] = nnneighbors[unique]

    return sites


def matching_sites(position, comparators, tol=1e-8):
    """ Get the indices of all points in a comparator list that are
    equal to a fractional coordinate (with a tolerance), taking into
    account periodic boundary conditions (adaptation from Pymatgen).

    Parameters:
      position: list (3,)
        Fractional coordinate to compare to list.
      comparators: list (3, n)
        Fractional coordinates to compare against.
      tol: float
        Absolute tolerance.

    Returns: list (n,)
        Indices of matches.
    """

    if len(comparators) == 0:
        return []

    fdist = comparators - position
    fdist -= np.round(fdist)
    return np.where((np.abs(fdist) < tol).all(axis=1))[0]


def ext_gcd(a, b):
    if b == 0:
        return 1, 0
    elif a % b == 0:
        return 0, 1
    else:
        x, y = ext_gcd(b, a % b)
        return y, x - y * (a // b)
