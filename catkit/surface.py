from . import utils
import numpy as np
from numpy.linalg import norm, pinv, solve
from ase.neighborlist import NeighborList
import networkx as nx
import networkx.algorithms.isomorphism as iso
from ase.build import rotate
from ase.constraints import FixAtoms
from scipy.spatial import Delaunay
from scipy.linalg import circulant
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

        Args:
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

        Args:
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

        # Orthogonolize the z-coordinate
        # Warning: bulk symmetry is lost at this point
        a1, a2, a3 = slab.cell
        a3 = (np.cross(a1, a2) * np.dot(a3, np.cross(a1, a2)) /
              norm(np.cross(a1, a2)) ** 2)
        slab.cell[2] = a3

        if self.vacuum:
            # Requires vacuum
            slab.center(vacuum=self.vacuum, axis=2)

            if primitive:
                slab = utils.get_primitive_cell(slab)
                zlayers = utils.get_unique_coordinates(
                    slab,
                    direct=False,
                    tag=True,
                    tol=self.tol
                )
                slab.rotate(slab.cell[0], 'x', rotate_cell=True)

                # spglib ocassionally returns a bimodal slab
                zpos = slab.get_scaled_positions().T[2]
                if zpos.max() > 0.9:
                    translate = slab.positions.T[2][zpos > 0.5].min()
                    slab.positions -= [0, 0, translate + self.tol]
                    slab.wrap(pbc=True)
                    slab.center(vacuum=self.vacuum, axis=2)

        elif primitive:
            raise(
                NotImplementedError,
                'Primitive slab generation requires vacuum'
            )

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
        """ Find the undercoordinated atoms at the upper and lower
        fraction of a given unit cell based on the bulk structure it
        origionated from.

        Assumes the xy-plane is perpendicualr to the miller index.

        Args:
          slab: ASE atoms-object
            The slab to find top layer atoms from.

        Returns: ndarray (n,), ndarray (m,)
          Array of atom indecies corresponding to the top and bottom
          layers of the slab.
        """

        if slab is None:
            slab = self.get_slab(primitive=True)

        ind, N = utils.get_voronoi_neighbors(self.bulk)

        radii = [self.bulk.get_distance(u, v, mic=True) for u, v in N.keys()]

        ind0, N0 = utils.get_cutoff_neighbors(slab, cutoff=max(radii))

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
    ):
        """ Helper function for getting the adsorption sites of a
        slab as defined by surface atom symmetries.

        Args:
          slab: ASE atoms-object
            The slab to find adsorption sites for.
          surface_sites: ndarray (n,)
            surface sites of the provided slab.
          unit_cell_rep: int
            Number of times to repeat the unit cell.
          tol: float
            Absolute tolerance for floating point errors.

        Returns: dict of 2 lists
          Dictionary of top, bridge, hollow, and 4fold sites containing
          two lists (positions and points) of the supercell.
        """

        if slab is None:
            slab = self.get_slab(primitive=True)

        if surface_sites is None:
            if self.surface_atoms is None:
                surface_sites, _ = self.get_surface_atoms(slab)

        sites = find_adsorption_sites(
            slab=slab,
            surface_sites=surface_sites,
            unit_cell_rep=5,
            tol=self.tol,
        )

        sites = get_reduced_sites(
            sites=sites,
            slab=slab,
            trim=[0., 4.],
            tol=self.tol,
        )

        if symmetry_reduced:
            sites = get_symmetric_sites(
                sites=sites,
                slab=slab,
                tol=self.tol,
            )

        return sites


def find_adsorption_sites(
        slab,
        surface_sites,
        unit_cell_rep=5,
        tol=1e-5,
):
    """ Find all bridge and hollow sites (3-fold and 4-fold) given an
    input slab based Delaunay triangulaton of surface atoms of a
    supercell.

    Args:
      slab: ASE atoms-object
        The slab to find adsorption sites for.
      surface_sites: ndarray (n,)
        surface sites of the provided slab.
      unit_cell_rep: int
        Number of times to repeat the unit cell.
      tol: float
        Absolute tolerance for floating point errors.

    Returns: dict of 2 lists
      Dictionary of top, bridge, hollow, and 4fold sites containing
      two lists (positions and points) of the supercell.
    """

    # Top sites projected into expanded unit cell
    n = len(slab)
    x = unit_cell_rep ** 2
    rtop = (np.arange(0, x*n, n).reshape(x, 1) + surface_sites).flatten()
    ratoms = slab * (unit_cell_rep, unit_cell_rep, 1)

    # Dict of 2 lists (positions, site_id)
    sites = {
        'top': [ratoms.positions[rtop], rtop],
        'bridge': [[], []],
        'hollow': [[], []],
        '4fold': [[], []],
    }

    dt = Delaunay(sites['top'][0].T[:2].T)
    neighbors = dt.neighbors
    simplices = dt.simplices

    bridge_neighbors = []

    for i, v in enumerate(simplices):

        cir = circulant(v)
        cornors = cir[2]
        edges = cir[:2]

        # Inner angle of each triangle corner
        vec = sites['top'][0][edges] - sites['top'][0][cornors]
        uvec = vec.T / norm(vec, axis=2).T
        angles = np.sum(uvec.T[0] * uvec.T[1], axis=1)

        # Angle types
        right = np.isclose(angles, 0)
        obtuse = (angles < -tol)

        rh_corner = cornors[right]
        edge_neighbors = neighbors[i][::-1]

        if obtuse.any():
            # Assumption: All simplices with obtuse angles
            # are irrelevent boundrys.
            continue

        bridge = np.sum(sites['top'][0][edges], axis=0) / 2.0

        # Looping through cornors allows for elimination of
        # redundent points, identification of 4-fold hollows,
        # and collection of bridge neigbors.
        for j, c in enumerate(cornors):

            edge = sorted(edges.T[j])

            if edge in sites['bridge'][1]:
                continue

            # Get the bridge neighbors (for adsorption vector)
            neighbor_simplex = simplices[edge_neighbors[j]]
            oc = list(set(neighbor_simplex) - set(edge))[0]

            # Right angles potentially indicate 4-fold hollow
            potential_hollow = sorted(edge + [c, oc])
            if c in rh_corner:

                if potential_hollow in sites['4fold'][1]:
                    continue

                # Assumption: If not 4-fold, this suggests
                # no hollow OR brige side is preasent.

                ovec = sites['top'][0][edge] - sites['top'][0][oc]
                ouvec = ovec / norm(ovec)
                oangle = np.dot(*ouvec)
                oright = np.isclose(oangle, 0)
                if oright:
                    sites['4fold'][0] += [bridge[j]]
                    sites['4fold'][1] += [rtop[potential_hollow].tolist()]

            else:
                sites['bridge'][0] += [bridge[j]]
                sites['bridge'][1] += [rtop[edge].tolist()]
                bridge_neighbors += [[c, oc]]

        if not right.any():
            hollow = np.average(sites['top'][0][v], axis=0)

            sites['hollow'][0] += [hollow]
            sites['hollow'][1] += [rtop[v].tolist()]

    # Set the reduced top sites
    sites['top'] = [slab.positions[surface_sites], surface_sites]

    # Convert lists to arrays and sort
    for k, v in sites.items():
        positions, points = v

        if len(positions) == 0:
            continue

        positions = np.array(positions)
        srt = np.lexsort((positions[:, 1], positions[:, 0]))

        sites[k][0] = positions[srt]
        sites[k][1] = np.array(points)[srt]

    return sites


def get_reduced_sites(
        sites,
        slab,
        trim=None,
        tol=1e-5,
):
    """ Reduce overlapping points via fractional coords. Intended
    for use after finding supercell sites.

    Args:
      sites: dict
        Dictionary of cartesian coordinates to provide reduced sites from.
        Must have the form {'site': [[positions], [points]]}.
      slab: ASE atoms-object
        slab to determine symmetry operations from.
      trim: list or tuple (2,)
        Lower and upper bounds of x and y coordinates to ommit.
      tol: float
        Absolute tolerance for floating point errors.

    Returns: dict of 2 lists
      Dictionary sites containing two lists (positions and points)
      of the supercell.
    """

    for k, v in sites.items():
        positions, points = v

        if len(positions) == 0:
            continue

        frac_coords = np.dot(positions, pinv(slab.cell))

        if isinstance(trim, (list, tuple)) and k != 'top':
            screen = (
                frac_coords.T[:2] > trim[0]).all(axis=0) & \
                (frac_coords.T[:2] < trim[1]).all(axis=0)

            frac_coords = frac_coords[screen]

        non_unique, unique_points, unique_positions = [], [], []
        for i, xyz in enumerate(frac_coords):
            if i in non_unique:
                continue

            non_unique += matching_sites(xyz, frac_coords).tolist()
            unique_positions += [xyz]
            unique_points += [points[i]]

        shift = [tol, tol, 0]
        unique_positions = np.array(unique_positions)
        unique_positions += shift

        unique_positions %= 1
        unique_positions -= shift

        sites[k][0] = np.dot(unique_positions, slab.cell)
        sites[k][1] = np.array(unique_points)

    return sites


def get_symmetric_sites(
        sites,
        slab,
        tol=1e-5
):
    """ Determine the symmetrically unique adsorption sites
    from a dictionary of possible sites for a given slab.

    Args:
      sites: dict
        Dictionary of cartesian coordinates to provide reduced sites from.
        Must have the form {'site': [[positions], [points]]}.
      slab: ASE atoms-object
        slab to determine symmetry operations from.
      tol: float
        Absolute tolerance for floating point errors.

    Returns: dict of 2 lists
      Dictionary sites containing two lists (positions and points)
      of the supercell.
    """

    symmetry = utils.get_symmetry(slab, tol=tol)
    rotations = symmetry['rotations']
    translations = symmetry['translations']

    for k, v in sites.items():
        positions, points = v

        if len(positions) == 0:
            continue

        frac_coords = np.dot(positions, pinv(slab.cell))

        # Convert to fractional
        unique_positions, unique_points = [], []
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
                unique_points += [points[i]]

        sites[k][0] = np.dot(unique_positions, slab.cell)
        sites[k][1] = np.array(unique_points)

    return sites


def matching_sites(position, comparators, tol=1e-8):
    """ Get the indices of all points in a comparator list that are
    equal to a fractional coordinate (with a tolerance), taking into
    account periodic boundary conditions. (adaptation from pymatgen).

    Args:
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
