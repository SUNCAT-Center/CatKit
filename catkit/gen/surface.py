from __future__ import division
from catkit import Gratoms
from . import utils
from . import adsorption
from . import geometry
import numpy as np
from numpy.linalg import norm, solve
from ase.build import rotate
from ase.constraints import FixAtoms
from itertools import product
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
                 miller_index,
                 layers,
                 vacuum=None,
                 fixed=None,
                 fix_stoichiometry=False,
                 attach_graph=True,
                 standardize_bulk=True,
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
        fix_stoichiometry : bool
            Constraints any slab generated to have the same
            stoichiometric ratio as the provided bulk.
        attach_graph : bool
            Attach the connectivity graph generated from the bulk structure.
            This is only necessary for fingerprinting and setting it to False
            can save time. Surface atoms will be found regardless.
        standardize_bulk : bool
            Covert the bulk input to its standard form before and
            produce the cleave from it. This is highly recommended as
            Miller indices are not defined for non-standard cells.
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
        self.attach_graph = attach_graph
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
            miller_index = np.round(miller_index)
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

        if not self.attach_graph:
            del new_bulk._graph

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
        ibasis : Gratoms object
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
        exbasis = ibasis * (1, 1, 2)

        connectivity = utils.get_voronoi_neighbors(exbasis)

        n = len(ibasis)
        diff = connectivity[:n, n:].sum(axis=1)
        surf_atoms = diff != 0

        if np.all(diff):
            warnings.warn(
                ("Your slab has no bulk atoms and may be too thin "
                 "to identify surface atoms correctly. This may cause "
                 "surface adsorption site identification to fail."))

        # TODO: Graph generation needs to go here once handling of
        # unit cell repetitions is implemented.
        scaled_zpositions = ibasis.get_scaled_positions()[:, 2] + self.tol
        scaled_zpositions = np.round(scaled_zpositions % 1 + self.tol, 4)

        indices = np.argwhere(surf_atoms).flatten()
        zcoords = scaled_zpositions - np.mean(scaled_zpositions)
        top = indices[zcoords[indices] >= 0]
        bottom = indices[zcoords[indices] < 0]

        ibasis.set_surface_atoms(top=top, bottom=bottom)

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
        slab = self.get_slab_basis(iterm).copy()

        # Trim the bottom of the cell, bulk symmetry will be lost
        if not self.fix_stoichiometry:
            zlayers = utils.get_unique_coordinates(slab, tol=self.tol)
            reverse_sort = np.sort(zlayers)[::-1]
            ncut = reverse_sort[:self.layers][-1] * slab.cell[2][2]

            zpos = slab.positions[:, 2]
            index = np.arange(len(slab))
            del slab[index[zpos - ncut < -self.tol]]

            slab.cell[2][2] -= ncut
            slab.translate([0, 0, -ncut])

        slab = self.set_size(slab, size)

        # Orthogonalize the z-coordinate
        # Breaks bulk periodicity in the c-basis
        slab.cell[2] = [0, 0, slab.cell[2][2]]
        slab.set_pbc([1, 1, 0])

        tl = np.min(slab.get_surface_atoms())
        translation = slab[tl].position.copy()
        translation[2] = 0
        slab.translate(-translation)
        slab.wrap()

        roundoff = np.isclose(slab.cell, 0)
        slab.cell[roundoff] = 0

        ind = np.lexsort(
            (slab.positions[:, 0],
             slab.positions[:, 1],
             slab.positions[:, 2]))
        slab = slab[ind]

        if self.fixed:
            tags = slab.get_tags()
            constraints = FixAtoms(mask=tags > (tags.max() - self.fixed))
            slab.set_constraint(constraints)

        if self.vacuum:
            slab.center(vacuum=self.vacuum, axis=2)

        self.slab = slab

        return slab

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
        output = adsorption.get_adsorption_sites(
            slab=slab, **kwargs)

        return output

    def set_size(self, slab, size):
        """Set the size of a slab based one of three methods.

        1. An integer value performs a search of valid matrix operations
        to perform on the ab-basis vectors to return a set which with
        a minimal sum of distances and an angle closes to 90 degrees.

        2. A list-like of length 2 will multiply the existing basis
        vectors by that amount.

        3. A list of two list-likes of length 2 will be interpreted
        as matrix notation to multiply the ab-basis vectors by.

        Parameters
        ----------
        slab : Atoms object
            Slab to be made into the requested size.
        size : int, list (2,) or (2, 2)
            Size of the unit cell to create as described above.

        Returns:
        --------
        supercell : Gratoms object
            Supercell of the requested size.
        """
        supercell = slab

        if isinstance(size, int):
            a = max(int(size / 2), 1) + size % 2
            T = np.mgrid[-a:a + 1, -a:a + 1].reshape(2, -1).T

            metrics = []
            for i, M in enumerate(product(T, repeat=2)):
                M = np.array(M)
                if ~np.isclose(abs(np.linalg.det(M)), size):
                    continue

                vector = np.dot(M.T, slab.cell[:2, :2])

                d = np.linalg.norm(vector, axis=1)
                angle = np.dot(vector[0], vector[1]) / np.prod(d)
                diff = np.diff(d)[0]

                # obtuse angle
                if angle < 0 or diff < 0:
                    continue

                metrics += [[d.sum(), angle, M]]

            if metrics:
                matrix = sorted(metrics, key=lambda x: (x[0], x[1]))[0][-1]
                supercell = transform_ab(supercell, matrix)

        elif isinstance(size, (list, tuple, np.ndarray)):
            size = np.array(size, dtype=int)

            if size.shape == (2,):
                supercell *= (size[0], size[1], 1)
            elif size.shape == (2, 2):
                supercell = transform_ab(supercell, size)

        if self.attach_graph:
            # TODO: Creating the graph at this point is not ideal.
            # Need to be able to handle expansion of the unit cell
            # before this can be moved back to basis creation
            n = len(supercell)
            exsupercell = supercell * (1, 1, 2)
            connectivity = utils.get_voronoi_neighbors(exsupercell)
            edges = utils.connectivity_to_edges(connectivity[:n, :n])
            supercell.graph.add_weighted_edges_from(edges, weight='bonds')

        return supercell


def transform_ab(slab, matrix, tol=1e-5):
    """Transform the slab basis vectors parallel to the z-plane
    by matrix notation. This can result in changing the slabs
    cell size. This can also result in very unusual slab dimensions,
    use with caution.

    Parameters
    ----------
    slab : Atoms object
        The slab to be transformed.
    matrix : array_like (2, 2)
        The matrix notation transformation of the a and b
        basis vectors.
    tol : float
        Float point precision tolerance.

    Returns
    -------
    slab : Atoms object
        Slab after transformation.
    """
    M = np.eye(3)
    M[:2, :2] = np.array(matrix).T
    newcell = np.dot(M, slab.cell)
    u, v = newcell[:2]
    nu, nv = norm(u), norm(v)

    # Ensure the longest axis is x without z-rotation.
    if nv > nu:
        matrix[:, [0, 1]] = matrix[:, [1, 0]]
        u *= -1

    if u[0] < 0 and v[1] > 0:
        matrix[[0, 1]] = matrix[[1, 0]]

    M[:2, :2] = np.array(matrix).T
    newcell = np.dot(M, slab.cell)

    scorners_newcell = np.array([
        [0, 0], [0, 0],
        [0, 1], [1, 1]])

    corners = np.dot(scorners_newcell, newcell[:2, :2])
    scorners = np.linalg.solve(slab.cell[:2, :2].T, corners.T).T

    rep = np.ceil(scorners.ptp(axis=0)).astype(int)

    slab *= (rep[0], rep[1], 1)
    slab.set_cell(newcell)

    coords = slab.get_scaled_positions()
    original_index = np.arange(coords.shape[0])
    periodic_match = original_index.copy()
    for i, j in enumerate(periodic_match):
        if i != j:
            continue

        matched = geometry.matching_sites(coords[i], coords)
        periodic_match[matched] = i

    repeated = np.where(periodic_match != original_index)
    del slab[repeated]

    # Align the longest basis vector with x
    slab.rotate(slab.cell[0], 'x', rotate_cell=True)

    if slab.cell[1][1] < 0:
        slab.cell[1] *= -1

    return slab


def ext_gcd(a, b):
    """Extension of greatest common divisor."""
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
