from __future__ import division
from catkit import Gratoms
from . import utils
from . import adsorption
import numpy as np
from numpy.linalg import norm, solve
from ase.build import rotate
from ase.constraints import FixAtoms
import networkx as nx
import warnings
import scipy
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
                 vacuum=None,
                 fixed=None,
                 standardize_bulk=True,
                 fix_stoichiometry=False,
                 tol=1e-8):
        """Generate a slab from a bulk atoms object.

        Return the miller indices associated with the users requested
        values. Follows the following steps:

        - Convert Miller-Bravais notation into standard Miller index.
        - Ensure the bulk cell is in its standard form.
        - Convert the indices to the cell for the primitive lattice.
        - Reduce the indices by their greatest common divisor.

        Parameters
        ----------
        bulk : Atoms object
            Bulk system to be converted into slab.
        miller_index : list (3,) or (4,)
            Miller index to construct surface from. If length 4, Miller-Bravais
            notation is assumed.
        layers : int
            Number of layers to include in the slab. A slab layer is defined
            as a unique z-coordinate.
        vacuum : float
            Angstroms of vacuum to apply to the slab.
        fixed : int
            Number of slab layers to constrain.
        standardize_bulk : bool
            Covert the bulk input to its standard form before and
            produce the cleave from it. This is highly recommended as
            Miller indices are not defined for non-standard cells.
        fix_stoichiometry : bool
            Constraints any slab generated to have the same
            stoichiometric ratio as the provided bulk.
        tol : float
            Tolerance for floating point rounding errors.
        """
        self.layers = layers
        self.vacuum = vacuum
        self.fixed = fixed
        self.tol = tol
        self.fix_stoichiometry = fix_stoichiometry
        self.unique_terminations = None
        self.standardized = standardize_bulk
        self.slab_basis = None
        self.slab = None

        miller_index = np.array(miller_index)
        if len(miller_index) == 4:
            miller_index[[0, 1]] -= miller_index[2]
            miller_index = np.delete(miller_index, 2)
        miller_index = (miller_index / list_gcd(miller_index)).astype(int)

        # Store the Miller indices associated with the standard cell.
        self.miller_index = miller_index

        if standardize_bulk:
            bulk = utils.get_spglib_cell(bulk)
            norm = get_normal_vectors(bulk)

            bulk = utils.get_spglib_cell(bulk, primitive=True)
            pnorm = get_normal_vectors(bulk)

            miller_index = np.dot(
                miller_index, np.dot(norm, np.linalg.inv(pnorm)))
            miller_index = (miller_index / list_gcd(miller_index)).astype(int)

            self._bulk = self.align_crystal(bulk, miller_index)

            z_planes = utils.get_unique_coordinates(self._bulk, tol=self.tol)
            self.unique_terminations = z_planes
            self.slab_basis = [None] * len(z_planes)
        else:
            self._bulk = self.align_crystal(bulk, miller_index)

    def align_crystal(self, bulk, miller_index):
        """Return a standardized unit cell from bulk unit cell. This
        standardization rotates the a and b basis vectors to be parallel
        to the Miller index.

        Parameters
        ----------
        bulk : Atoms object
            Bulk system to be standardized.
        miller_index : list (3,)
            Miller indices to align with the basis vectors.

        Returns
        -------
        new_bulk : Gratoms object
            Standardized bulk unit cell.
        """
        del bulk.constraints

        if len(np.nonzero(miller_index)[0]) == 1:
            mi = max(np.abs(miller_index))
            new_lattice = scipy.linalg.circulant(
                miller_index[::-1] / mi).astype(int)
        else:
            h, k, l = miller_index
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
            new_lattice = np.array([c1, c2, c3])

        scaled = solve(new_lattice.T, bulk.get_scaled_positions().T).T
        scaled -= np.floor(scaled + self.tol)

        new_bulk = Gratoms(
            positions=bulk.positions,
            numbers=bulk.get_atomic_numbers(),
            pbc=True)

        new_bulk.set_scaled_positions(scaled)
        new_bulk.set_cell(np.dot(new_lattice, bulk.cell), scale_atoms=True)

        # Align the longest of the ab basis vectors with x
        d = norm(new_bulk.cell[:2], axis=1)
        if d[1] > d[0]:
            new_bulk.cell[[0, 1]] = new_bulk.cell[[1, 0]]
        a = new_bulk.cell[0]
        a3 = np.cross(a, new_bulk.cell[1]) / np.max(d)
        rotate(new_bulk, a3, (0, 0, 1), a, (1, 0, 0))

        # Ensure the remaining basis vectors are positive in their
        # corresponding axis
        for i in range(1, 3):
            if new_bulk.cell[i][i] < 0:
                new_bulk.cell[i] *= -1

        return new_bulk

    def get_unique_terminations(self):
        """Determine the fractional coordinate shift that will result in
        a unique surface termination. This is not required if bulk
        standardization has been performed, since all available z shifts will
        result in a unique termination for a primitive cell.

        Returns
        -------
        unique_shift : array (n,)
            Fractional coordinate shifts which will result in unique
            terminations.
        """
        if self.unique_terminations is not None:
            return self.unique_terminations

        rotations, translations = utils.get_symmetry(self._bulk, tol=self.tol)

        # Find all symmetries which are rotations about the z-axis
        zsym = np.abs(rotations)
        zsym[:, 2, 2] -= 1
        zsym = zsym[:, [0, 1, 2, 2, 2], [2, 2, 2, 0, 1]]
        zsym = np.argwhere(zsym.sum(axis=1) == 0)

        itol = self.tol ** -1
        ztranslations = np.floor(translations[zsym, -1] * itol) / itol
        z_symmetry = np.unique(ztranslations)

        zcoords = utils.get_unique_coordinates(self._bulk, tol=self.tol)
        zdiff = np.cumsum(np.diff(zcoords))
        zdiff = np.floor(zdiff * itol) / itol
        unique_shift = np.argwhere(zdiff < z_symmetry[1]) + 1
        unique_shift = np.append(0, zcoords[unique_shift])

        self.unique_terminations = unique_shift
        self.slab_basis = [None] * len(unique_shift)

        return unique_shift

    def get_slab_basis(self, iterm):
        """Return a list of all terminations which have been properly shifted
        and with an appropriate number of layers added. This function is mainly
        for performance, to prevent looping over other operations which are not
        related the size of the slab.

        Only produces the terminations requested.

        Parameters
        ----------
        iterm : int
            Index of the slab termination to return.

        Returns
        -------
        ibasis :Gratoms object
            Prepared, ith basis.
        """
        terminations = self.get_unique_terminations()

        if self.slab_basis[iterm] is not None:
            ibasis = self.slab_basis[iterm].copy()
            return ibasis

        _basis = self._bulk.copy()
        if self.standardized:
            bulk_layers = self.unique_terminations
        else:
            bulk_layers = utils.get_unique_coordinates(_basis, tol=self.tol)

        zshift = terminations[iterm]
        ibasis = _basis.copy()
        ibasis.translate([0, 0, -zshift])
        ibasis.wrap(pbc=True)

        minimum_repetitions = np.ceil(self.layers / len(bulk_layers))
        ibasis *= (1, 1, int(minimum_repetitions))

        utils.get_unique_coordinates(
            ibasis, tag=True, tol=self.tol)

        self.slab_basis[iterm] = ibasis

        return ibasis

    def get_slab(self, size=(1, 1), iterm=0):
        """Generate a slab from the bulk structure. This function is meant
        specifically selection of an individual termination or enumeration
        through surfaces of various size.

        This function will orthogonalize the c basis vector and align it
        with the z-axis which breaks bulk symmetry along the z-axis.

        Parameters
        ----------
        size : tuple (2,)
            Repeat the x and y lattice vectors by the indicated
            dimensions
        iterm : int
            A termination index in reference to the list of possible
            terminations.

        Returns
        -------
        slab : Atoms object
            The modified basis slab produced based on the layer specifications
            given.
        """
        slab = self.get_slab_basis(iterm)

        # Orthogonalize the z-coordinate
        slab.cell[2] = [0, 0, slab.cell[2][2]]

        # Trim the bottom of the cell
        if not self.fix_stoichiometry:
            zlayers = utils.get_unique_coordinates(slab, tol=self.tol)
            reverse_sort = np.sort(zlayers)[::-1]
            ncut = reverse_sort[:self.layers][-1] * slab.cell[2][2]

            zpos = slab.positions[:, 2]
            index = np.arange(len(slab))
            del slab[index[zpos - ncut < -self.tol]]

            slab.cell[2][2] -= ncut
            slab.translate([0, 0, -ncut])

        slab *= (size[0], size[1], 1)

        translation = slab[0].position.copy()
        translation[2] = 0
        slab.translate(-translation)
        slab.wrap()

        ind = np.lexsort(
            (slab.positions[:, 0],
             slab.positions[:, 1],
             slab.positions[:, 2]))
        slab = slab[ind]
        slab.set_pbc([1, 1, 0])

        if self.fixed:
            tags = slab.get_tags()
            constraints = FixAtoms(mask=tags > (tags.max() - self.fixed))
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
        bulk_con = utils.get_voronoi_neighbors(self._bulk)
        d = self._bulk.get_all_distances(mic=True)
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
        bulk_con = utils.get_voronoi_neighbors(self._bulk)
        d = self._bulk.get_all_distances(mic=True)

        # This is potentially problematic for very small unit cells
        # as it does not take periodicity into account
        maxd = d[bulk_con > 0].max()
        bulk_con = utils.get_cutoff_neighbors(self._bulk, cutoff=maxd)

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
    if isinstance(values[0], float):
        values = np.array(values, dtype=int)

    gcd_func = np.frompyfunc(gcd, 2, 1)
    _gcd = np.ufunc.reduce(gcd_func, values)

    return _gcd


def get_normal_vectors(atoms):
    """Return the normal vectors to a atoms unit cell."""
    rotation1 = np.roll(atoms.cell, 1, axis=0)
    rotation2 = np.roll(atoms.cell, -1, axis=0)
    normal_vectors = np.cross(rotation1, rotation2, axis=0)

    return normal_vectors
