from __future__ import division
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import absolute_import
from builtins import int
from builtins import round
from builtins import range
from future import standard_library
standard_library.install_aliases()
from . import utils
from numpy.linalg import norm
import numpy as np
import math
from ase import Atoms
from ase.neighborlist import NeighborList
import networkx as nx
import networkx.algorithms.isomorphism as iso
from ase.constraints import FixAtoms


class SlabGenerator(object):
    """ Class for generation of slab unit cells from bulk unit cells.
    """

    def __init__(
            self,
            bulk,
            miller_index=[1, 1, 1],
            tol=1e-5):
        """ Generate a slab from an ASE bulk atoms-object.

        Parameters:
          bulk: ASE atoms-object
            Bulk structure to produce the slab from.

          miller_index: list (3,)
            Miller index to construct surface from.

          tol: float
            Tolerance for floating point rounding errors.
        """

        self.bulk = bulk
        self.miller_index = np.array(miller_index)
        self.tol = tol
        self.unique_terminations = None
        self.basis = None

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

        gcd = math.gcd(abs(self.miller_index[0]), abs(self.miller_index[1]))
        gcd = math.gcd(gcd, abs(self.miller_index[2]))

        self.miller_index = np.round(self.miller_index / gcd).astype(int)

        # Find surafce lattice vectors in int-math
        p, q = self._get_surface_vectors()

        # Get surface vectors in floating-math
        r = np.cross(p, q)
        p = np.dot(p, self.bulk.cell)
        q = np.dot(q, self.bulk.cell)
        r = np.dot(r, self.bulk.cell)

        slab = Atoms(cell=np.vstack([p, q, r]))

        # get the number of atoms in surface-basis
        natoms = int(round(
            slab.get_volume() / self.bulk.get_volume()))
        bulk_to_surface = np.dot(self.bulk.cell, np.linalg.inv(slab.cell))
        direct = self.bulk.get_scaled_positions()

        # get the required number of atoms
        for atom in self.bulk:

            P_surface = np.dot(direct[atom.index], bulk_to_surface)

            left = -self.tol - P_surface
            right = (1 - self.tol) - P_surface

            N = np.array([-1, 1])
            counter = 0

            # Loop until we get required number of atoms
            while counter != natoms:
                tmp_atoms = Atoms(cell=slab.get_cell())
                counter = 0
                N *= 2

                for t1 in range(N[0], N[1]):
                    for t2 in range(N[0], N[1]):
                        for t3 in range(N[0], N[1]):
                            temp = sum(bulk_to_surface[:, 0] * [t1, t2, t3])
                            if temp < left[0] or temp > right[0]:
                                continue

                            temp = sum(bulk_to_surface[:, 1] * [t1, t2, t3])
                            if temp < left[1] or temp > right[1]:
                                continue

                            temp = sum(bulk_to_surface[:, 2] * [t1, t2, t3])
                            if temp < left[2] or temp > right[2]:
                                continue

                            new_pos = np.dot(
                                np.add(direct[atom.index], [t1, t2, t3]),
                                bulk_to_surface)

                            slab_atom = Atoms(
                                atom.symbol,
                                cell=slab.get_cell())
                            slab_atom.set_scaled_positions(new_pos)
                            tmp_atoms += slab_atom

                            counter += 1

            slab += tmp_atoms

        # Rotate to align with x-axis with the surface  and z axis
        # perpendicular.
        a1, a2, a3 = slab.cell
        slab.set_cell([
            (norm(a1), 0, 0),
            (np.dot(a1, a2) / norm(a1),
             np.sqrt(norm(a2)**2 - (np.dot(a1, a2) / norm(a1))**2), 0),
            (0, 0, norm(a3))],
                      scale_atoms=True)

        slab.pbc = (True, True, False)

        self.basis = slab

        return slab

    def get_unique_terminations(
            self,
            primitive_cell=True):
        """ Return smallest unit cell corresponding to given surface and
        unique surface terminations based on symmetry and nearest neighbors.

        This function is mostly functionality better suited for a bulk
        enumerator class.

        Parameters:
          primitive_cell: bool
            Attempt to convert the unit cell into a primitive cell. This is
            temporary and should be moved to bulk enumeration functionality.

        Returns:
          unique_terminations: list
            Unique terminations of a surface.
        """

        # Get surface basis unit cell
        if self.basis is None:
            self.build_basis()

        # FUTURE: This should be an bulk generator operation.
        # For now, this redefines the internal miller_index
        if primitive_cell:
            b_base = [
                np.cross(self.bulk.cell[1], self.bulk.cell[2]),
                np.cross(self.bulk.cell[2], self.bulk.cell[0]),
                np.cross(self.bulk.cell[0], self.bulk.cell[1])]

            pbulk = utils.get_primitive_cell(
                self.bulk,
                tol=self.tol)

            b_prim = [
                np.cross(pbulk.cell[1], pbulk.cell[2]),
                np.cross(pbulk.cell[2], pbulk.cell[0]),
                np.cross(pbulk.cell[0], pbulk.cell[1])]

            miller_index = np.dot(
                self.miller_index,
                np.dot(b_base, np.linalg.inv(b_prim)))

            self.bulk = pbulk
            self.miller_index = np.round(miller_index).astype(int)

        unique_terminations = self._get_unique_terminations()
        self.unique_terminations = unique_terminations

        return unique_terminations

    def _get_unique_terminations(self):
        """ Return smallest unit cell corresponding to given surface and
        unique surface terminations based on symmetry and nearest neighbors.

        Returns:
          unique_terminations: list
            Unique terminations of a surface.
        """

        # Find all different planes as simply different z-coordinates
        z_planes = utils.get_unique_coordinates(self.basis, tol=self.tol)

        # now get the symmetries of lattice
        symmetry = utils.get_symmetry(self.basis, tol=self.tol)
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
            tmp_slab = self.basis.copy()
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

        return unique_terminations

    def _get_surface_vectors(self):
        """ Find integral inplane surface-vectors correspong to given integral plane.
        Given integral vectors, find p and q of minimum length such that p * q
        remains unchanged

        This function is meant for internal use.

        Returns:
          vectors: ndarray (2,)
            Set of integral-triplets representing the  minimum length
            surface vectors.
        """

        limits = np.array([-1, 1])
        while True:
            vectors = utils.vector_search(self.miller_index, limits)

            if vectors is not None:
                break

            limits *= 2

        # Find vectors with minimum lengths
        p, q = vectors
        len_p = sum(p ** 2)
        len_q = sum(q ** 2)

        flipped = False
        if len_p > len_q:
            p, q = q, p
            flipped = True

        alpha = int(round(np.dot(p, q) / np.dot(p, p)))
        q = q - alpha * p

        # flip back so that cross product is one
        if flipped:
            p, q = q, p

        return np.array([p, q])

    def get_slab(
            self,
            layers=3,
            fixed=2,
            iterm=0):
        """ Handler function for generating a slab object with a certain number
         of layers.

        Parameters:
          layers: int
            Number of layers to include in the slab.

          fixed: int
            Number of layers to fix in the slab.

          iterm: int
            A termination index in reference to the list of possible
            terminations.

        Returns:
          slab: ASE atoms-object
            The modified basis slab produced based on the layer specifications
            given.
        """

        if self.basis is None:
            slab = self.build_basis()
        slab = self.basis.copy()

        if iterm:
            terminations = self.unique_terminations
            zshift = terminations[iterm]

            slab.translate([0, 0, -zshift])
            slab.wrap(pbc=True)

        slab = self.set_layers(slab, layers, self.tol)
        slab = self.set_fixed_layers(slab, fixed, self.tol)

        return slab

    def set_layers(
            self,
            slab,
            layers,
            atol=1e-5):
        """ Convert a basis slab into the requested number of layers. Each
        layer is identified as a unique z-coordinate.

        This needs to be an operation on a slab class at some point.

        Parameters:
          slab: ASE atom-object
            A basis slab with an unspecified number of layers to
            modify.

          layers: int
            Number of layers to include in the slab.

          atol: float
            Tolerance for dealing with floating point rounding errors.

        Returns:
          slab: ASE atoms-object
            A modified slab with the desired number of layers.
        """

        zlayers = utils.get_unique_coordinates(slab, tol=atol)

        z_repetitions = math.ceil(layers / len(zlayers))
        slab *= (1, 1, z_repetitions)

        zlayers = utils.get_unique_coordinates(
            slab,
            direct=False,
            tol=atol)
        ncut = sorted(zlayers)[::-1][:layers][-1]

        zpos = slab.positions[:, 2]
        index = np.arange(len(slab))
        del slab[index[zpos - ncut < -atol]]

        slab.cell[2][2] -= ncut
        slab.translate([0, 0, -ncut])

        return slab

    def set_fixed_layers(
            self,
            slab,
            fixed,
            atol=1e-5):
        """ Convert a basis slab into the requested number of fixed layers.
        Each layer is identified as a unique z-coordinate.

        This needs to be an operation on a slab class at some point.

        Parameters:
          slab: ASE atom-object
            A basis slab with an unspecified number of layers to
            modify.

          fixed: int
            Number of layers to fix in the slab.

          atol: float
            Tolerance for dealing with floating point rounding errors.

        Returns:
          slab: ASE atoms-object
            A modified slab with the desired number of fixed layers.
        """

        del slab.constraints

        zlayers = utils.get_unique_coordinates(
            slab,
            direct=False,
            tol=atol)
        ncut = sorted(zlayers)[fixed]

        zpos = slab.positions[:, 2]
        index = np.arange(len(slab))

        constraints = FixAtoms(indices=index[zpos - ncut < -atol])
        slab.set_constraint(constraints)

        return slab
