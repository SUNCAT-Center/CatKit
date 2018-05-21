from __future__ import division
from catkit import Gratoms
from . import utils
from . import adsorption
from . import geometry
import numpy as np
from numpy.linalg import norm, solve
from ase.build import rotate
from ase.constraints import FixAtoms
import networkx as nx
import warnings
from scipy.linalg import circulant
try:
    from math import gcd
except ImportError:
    from fractions import gcd


class SlabGenerator(object):
    """Class for generation of slab unit cells from bulk unit cells.

    Many surface operations rely upon / are made easier through the
    bulk basis cell they are created from. The SlabGenerator class
    is designed to house these operations.
    """

    def __init__(self,
                 bulk,
                 miller_index=(1, 1, 1),
                 layers=4,
                 min_width=None,
                 fixed=0,
                 vacuum=0,
                 fix_stoichiometry=False,
                 tol=1e-8):
        """Generate a slab from a bulk atoms object.

        Parameters
        ----------
        bulk : atoms object
            Bulk structure to produce the slab from.
        miller_index : list (3,) or (4,)
            Miller index to construct surface from. If length is 4, then
            Miller-Bravis indices are assumed.
        layers : int
            Number of layers to include in the slab.
        min_width : float
            Minimum width of slabs produce (Ang). Will override layers.
        fixed : int
            Number of layers to fix in the slab.
        vacuum : float
            Angstroms of vacuum to add to the slab.
        fix_stoichiometry : bool
            Constraints any slab generated to have the same
            stoichiometry as the provided bulk.
        tol : float
            Tolerance for floating point rounding errors.
        """
        miller_index = np.array(miller_index)

        if len(miller_index) == 4:
            miller_index[[0, 1]] -= miller_index[2]
            miller_index = np.delete(miller_index, 2)
        miller_index = (miller_index / list_gcd(miller_index)).astype(int)

        self.miller_index = miller_index
        self.layers = layers
        self.fixed = fixed
        self.min_width = min_width
        self.vacuum = vacuum
        self.tol = tol

        self.fix_stoichiometry = fix_stoichiometry
        self.unique_terminations = None
        self.slab = None

        self._basis = self._build_basis(bulk)

    def _build_basis(self, bulk):
        """Get the basis unit cell from bulk unit cell. This
        basis is effectively the same as the bulk, but rotated such
        that the z-axis is aligned with the surface termination.

        The basis is stored separately from the slab generated and is
        only intended for internal use.

        Returns
        -------
        basis : atoms object
            The basis slab corresponding to the provided bulk.
        """
        del bulk.constraints

        if len(np.nonzero(self.miller_index)[0]) == 1:
            mi = max(np.abs(self.miller_index))
            basis = circulant(self.miller_index[::-1] / mi).astype(int)
        else:
            h, k, l = self.miller_index
            p, q = ext_gcd(k, l)
            a1, a2, a3 = bulk.cell

            k1 = np.dot(p * (k * a1 - h * a2) + q * (l * a1 - h * a3),
                        l * a2 - k * a3)
            k2 = np.dot(l * (k * a1 - h * a2) - k * (l * a1 - h * a3),
                        l * a2 - k * a3)

            if abs(k2) > self.tol:
                i = -int(np.round(k1 / k2))
                p, q = p + i * l, q - i * k

            a, b = ext_gcd(p * k + q * l, h)

            c1 = (p * k + q * l, -p * h, -q * h)
            c2 = np.array((0, l, -k)) // abs(gcd(l, k))
            c3 = (b, a * p, a * q)

            basis = np.array([c1, c2, c3])

        basis_atoms = Gratoms(
            positions=bulk.positions,
            numbers=bulk.get_atomic_numbers(),
            cell=bulk.cell,
            pbc=True)

        scaled = solve(basis.T, basis_atoms.get_scaled_positions().T).T
        scaled -= np.floor(scaled + self.tol)
        basis_atoms.set_scaled_positions(scaled)
        basis_atoms.set_cell(np.dot(basis, basis_atoms.cell), scale_atoms=True)

        a1, a2, a3 = basis_atoms.cell
        n1 = np.cross(a1, a2)
        a3 = n1 / norm(n1)
        rotate(basis_atoms, a3, (0, 0, 1), a1, (1, 0, 0))

        return basis_atoms

    def get_unique_terminations(self):
        """Return smallest unit cell corresponding to given surface and
        unique surface terminations based on symmetry and nearest neighbors.

        Returns
        -------
        unique_shift : list
            Unique terminations of a surface.
        """
        z_planes = utils.get_unique_coordinates(self._basis, tol=self.tol)

        symmetry = utils.get_symmetry(self._basis, tol=self.tol)
        rotations = symmetry['rotations']
        translations = symmetry['translations']

        # Find all symmetries which are rotations about the z-axis
        z_symmetry = []
        for i, rotation in enumerate(rotations):
            if (abs(rotation[2][0]) < self.tol
                    and abs(rotation[2][1]) < self.tol
                    and abs(rotation[0][2]) < self.tol
                    and abs(rotation[1][2]) < self.tol
                    and abs(rotation[2][2] - 1.0) < self.tol):

                if not np.isclose(
                        translations[i][2], z_symmetry, rtol=self.tol).any():
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

        self.unique_terminations = unique_shift

        return unique_shift

    def get_slab(self, size=(1, 1), root=None, iterm=None, primitive=False):
        """Generate a slab object with a certain number of layers.

        Parameters
        ----------
        size : tuple (2,)
            Repeat the x and y lattice vectors by the indicated
            dimensions
        root : int
            Produce a slab with a primitive a1 basis vector multiplied
            by the square root of a provided value. Uses primitive unit cell.
        iterm : int
            A termination index in reference to the list of possible
            terminations.
        primitive : bool
            Whether to reduce the unit cell to its primitive form.

        Returns
        -------
        slab : atoms object
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
            slab, direct=False, tol=self.tol)

        if self.min_width:
            width = slab.cell[2][2]
            z_repetitions = np.ceil(width / len(zlayers) * self.min_width)
        else:
            z_repetitions = np.ceil(self.layers / len(zlayers))

        slab *= (1, 1, int(z_repetitions))

        if primitive or root:
            if self.vacuum:
                slab.center(vacuum=self.vacuum, axis=2)
            else:
                raise ValueError('Primitive slab generation requires vacuum')

            nslab = utils.get_primitive_cell(slab)

            if nslab is not None:
                slab = nslab

                # spglib occasionally returns a split slab
                zpos = slab.get_scaled_positions()
                if zpos[:, 2].max() > 0.9 or zpos[:, 2].min() < 0.1:
                    zpos[:, 2] -= 0.5
                    zpos[:, 2] %= 1
                    slab.set_scaled_positions(zpos)
                    slab.center(vacuum=self.vacuum, axis=2)

                # For hcp(1, 1, 0), primitive alters z-axis
                d = norm(slab.cell, axis=0)
                maxd = np.argwhere(d == d.max())[0][0]
                if maxd != 2:
                    slab.rotate(slab.cell[maxd], 'z', rotate_cell=True)
                    slab.cell[[maxd, 2]] = slab.cell[[2, maxd]]
                    slab.cell[maxd] = -slab.cell[maxd]
                    slab.wrap(pbc=True)

                slab.rotate(slab.cell[0], 'x', rotate_cell=True)

        # Orthogonalize the z-coordinate
        # Warning: bulk symmetry is lost at this point
        a, b, c = slab.cell
        nab = np.cross(a, b)
        c = (nab * np.dot(c, nab) / norm(nab)**2)
        slab.cell[2] = c

        # Align the longest remaining basis vector with x
        vdist = norm(slab.cell[:2], axis=1)
        if vdist[1] > vdist[0]:
            slab.rotate(slab.cell[1], 'x', rotate_cell=True)
            slab.cell[0] *= -1
            slab.cell[[0, 1]] = slab.cell[[1, 0]]

        # Enforce that the angle between basis vectors is acute.
        if slab.cell[1][0] < 0:
            slab.rotate(slab.cell[2], '-z')
            slab.cell *= [[1, 0, 0], [-1, 1, 0], [0, 0, 1]]

        if root:
            roots, vectors = root_surface_analysis(
                slab, return_vectors=True)

            if root not in roots:
                raise ValueError(
                    'Requested root structure unavailable for this system.'
                    'Try: {}'.format(roots))

            vect = vectors[np.where(root == roots)][0]
            slab = self.root_surface(slab, root, vect)

        # Get the direct z-coordinate of the requested layer
        zlayers = utils.get_unique_coordinates(
            slab, direct=False, tag=True, tol=self.tol)

        if not self.fix_stoichiometry:
            reverse_sort = np.sort(zlayers)[::-1]

            if self.min_width:
                n = np.where(zlayers < self.min_width, 1, 0).sum()
                ncut = reverse_sort[n]
            else:
                ncut = reverse_sort[:self.layers][-1]

            zpos = slab.positions[:, 2]
            index = np.arange(len(slab))
            del slab[index[zpos - ncut < -self.tol]]

            slab.cell[2][2] -= ncut
            slab.translate([0, 0, -ncut])

        slab *= (size[0], size[1], 1)
        tags = slab.get_tags()

        m = np.where(tags == 1)[0][0]
        translation = slab[m].position.copy()
        translation[2] = 0
        slab.translate(-translation)
        slab.wrap()

        ind = np.lexsort(
            (slab.positions[:, 0],
             slab.positions[:, 1],
             slab.positions[:, 2]))

        slab = Gratoms(
            positions=slab.positions[ind],
            numbers=slab.numbers[ind],
            cell=slab.cell,
            pbc=[1, 1, 0],
            tags=tags[ind]
        )

        fix = tags.max() - self.fixed
        constraints = FixAtoms(indices=[a.index for a in slab if a.tag > fix])
        slab.set_constraint(constraints)

        if self.vacuum:
            slab.center(vacuum=self.vacuum, axis=2)

        self.slab = slab

        return slab

    def get_graph_from_bulk(self, slab, attach=False):
        """Return the surface connectivity of a slab based on information
        from the bulk basis is was constructed from.

        Assumes the xy-plane is perpendicular to the miller index.

        Parameters
        ----------
        slab : atoms object
            Atoms to find graph connectivity for.

        Returns
        -------
        surf_con : ndarray (n, n)
            Connectivity matrix of the surface atoms with periodic boundary
            conditions.
        """
        bulk_con = utils.get_voronoi_neighbors(self._basis)
        d = self._basis.get_all_distances(mic=True)
        maxd = d[bulk_con > 0].max()

        surf_con = utils.get_cutoff_neighbors(slab, cutoff=maxd)

        if attach:
            G = nx.MultiGraph(surf_con)
            slab.graph.add_edges_from(G.edges(data=True))

        return surf_con

    def get_voronoi_surface_atoms(self, slab, attach_graph=True):
        """Find the under-coordinated atoms at the upper and lower
        fraction of a given unit cell based on the bulk structure it
        originated from.

        Assumes the xy-plane is perpendicular to the miller index.

        Parameters
        ----------
        slab : atoms object
            The slab to find top layer atoms from.
        attach_graph : bool
            Store the slabs graph information in the passed  slab object.

        Returns
        -------
        top : ndarray (n,)
            Array of atom indices corresponding to the top layer of the slab.
        bottom : ndarray (m,)
            Array of atom indices corresponding to the bottom layer of
            the slab.
        """
        bulk_con = utils.get_voronoi_neighbors(self._basis)
        d = self._basis.get_all_distances(mic=True)

        # This is potentially problematic for very small unit cells
        # as it does not take periodicity into account
        maxd = d[bulk_con > 0].max()
        bulk_con = utils.get_cutoff_neighbors(self._basis, cutoff=maxd)

        bulk_degree = bulk_con.sum(axis=0)
        surf_con = utils.get_cutoff_neighbors(slab, cutoff=maxd)
        surf_degree = surf_con.sum(axis=0)

        if attach_graph:
            edges = utils.connectivity_to_edges(surf_con)
            slab.graph.add_weighted_edges_from(edges, weight='bonds')

        # Expand bulk topology to match number of atoms in slab
        rep = np.ceil(len(surf_degree) / len(bulk_degree))
        bulk_degree = np.repeat(bulk_degree, rep)

        diff_degree = surf_degree - bulk_degree[:len(surf_degree)]
        surf_atoms = np.nonzero(diff_degree)[0]
        bulk_atoms = set(np.arange(len(slab))) - set(surf_atoms)

        if len(bulk_atoms) == 0:
            warnings.warn(
                ("Your slab has no bulk atoms and may be too thin "
                 "to identify surface atoms correctly. This may cause "
                 "surface adsorption site identification to fail."))

        hwp = slab.positions[surf_atoms] - slab.get_center_of_mass()
        top = surf_atoms[hwp.T[2] > 0]
        bottom = surf_atoms[hwp.T[2] < 0]

        return top, bottom

    def adsorption_sites(self, slab, **kwargs):
        """Helper function to return the adsorption sites of the provided slab.

        Parameters
        ----------
        slab : atoms object
            The slab to find adsorption sites for. Assumes you are using
            the same basis.

        Returns
        -------
        coordinates : ndarray (n,)
            Coordinates of the adsorption sites
        connectivity : ndarray (n,)
            Connectivity of the adsorption sites
        """
        if slab != self.slab or slab.get_surface_atoms() is None:
            surface_sites = self.get_voronoi_surface_atoms(slab)[0]
            slab.set_surface_atoms(surface_sites)
            self.slab = slab

        output = adsorption.get_adsorption_sites(
            slab=slab, **kwargs)

        return output

    def root_surface(self, slab, root, vector=None):
        """Creates a cell from a primitive cell that repeats along the x and y
        axis in a way consistent with the primitive cell, that has been cut
        to have a side length of *root*.

        *primitive cell* should be a primitive 2d cell of your slab, repeated
        as needed in the z direction.

        *root* should be determined using an analysis tool such as the
        root_surface_analysis function, or prior knowledge. It should always
        be a whole number as it represents the number of repetitions.
        """
        atoms = slab.copy()

        xynorm = norm(atoms.cell[0:2, 0:2], axis=1)
        vectors = atoms.cell[0:2, 0:2] / xynorm[0]
        mvectors = norm(vectors, axis=1)

        vnorm = norm(vector)
        angle = -np.arctan2(vector[1], vector[0])
        scale = vnorm / mvectors[0]

        maxn = int(np.ceil(scale)) * 2
        atoms *= (maxn, maxn, 1)
        atoms.cell[0:2, 0:2] = slab.cell[0:2, 0:2] * scale
        atoms.rotate((angle * 180) / np.pi, 'z')

        coords = atoms.get_scaled_positions()
        original_index = np.arange(coords.shape[0])
        periodic_match = original_index.copy()
        for i, j in enumerate(periodic_match):
            if i != j:
                continue

            matched = geometry.matching_sites(coords[i], coords)
            periodic_match[matched] = i

        repeated = np.where(periodic_match != original_index)
        del atoms[repeated]

        atoms.wrap()
        ind = np.lexsort(
            (atoms.positions[:, 0],
             atoms.positions[:, 1],
             atoms.positions[:, 2]))
        atoms = Gratoms(
            positions=atoms.positions[ind],
            numbers=atoms.numbers[ind],
            cell=atoms.cell,
            pbc=atoms.pbc,
            tags=atoms.get_tags()
        )

        assert(len(atoms) == len(slab) * root)

        return atoms


def root_surface_analysis(slab, max_cell=20, return_vectors=True):
    """A tool to analyze a slab and look for valid roots that exist, up to
    the given cell size. This is adapted from ASE.

    Parameters
    ----------
    slab : atoms object
        The slab to find adsorption sites for. Assumes you are using
        the same basis.
    max_cell : bool
        The maximum cell size to search for root structures within.
    return_vectors : bool
        Return the vectors associated with the valid root instances.

    Returns
    -------
    valid_roots : ndarray (n,)
        Return the integers which represent valid root distances.
        to create a unit cell from.
    bvectors :  ndarray (n, 2)
        The xy vectors associated with each of the valid integers.
    """
    xnorm = norm(slab.cell[0][0:2])
    vectors = slab.cell[0:2, 0:2] / xnorm

    mvectors = norm(vectors, axis=1)
    space = np.ceil(max_cell * 1.2 / mvectors)

    grid = np.mgrid[0:space[0], 0:space[1]]
    vect = np.dot(grid.T, vectors).reshape(-1, 2)[1:]

    dist = np.around(np.sum(vect ** 2, axis=1),
                     -int(np.log10(1e-8)))

    integers = np.equal(dist % 1, 0)
    select = np.where((dist <= max_cell) & integers)
    valid_roots, unique = np.unique(dist[select], return_index=True)
    valid_roots = valid_roots.astype(int)

    if return_vectors:
        bvectors = vect[select][unique]
        return valid_roots, bvectors

    return valid_roots


def ext_gcd(a, b):
    """Extension of greatest common denominator."""
    if b == 0:
        return 1, 0
    elif a % b == 0:
        return 0, 1
    else:
        x, y = ext_gcd(b, a % b)
        return y, x - y * (a // b)


def list_gcd(values):
    """Return the greatest common divisor of a list of values."""
    gcd_func = np.frompyfunc(gcd, 2, 1)
    _gcd = np.ufunc.reduce(gcd_func, values)

    return _gcd
