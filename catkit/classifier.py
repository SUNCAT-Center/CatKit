from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
import networkx as nx
from ase.neighborlist import NeighborList as NL
from ase.data import atomic_numbers as an
from ase.data import covalent_radii as r
from catkit.util import get_neighbors


class Classifier(object):
    """ Class for classification of various aspects of an an atomic
    unit cell.

    Currently, a tool for classification of adsorbates on surface
    environments and the active sites they rest on.
    """

    def __init__(
            self,
            atoms):
        """ Return unique coordinate values of a given atoms object
        for a specified axis.

        Parameters:
          atoms: ASE atoms-object
        """

        self.atoms = atoms
        self.neighbor_list = None
        self.molecules = None

    def _build_neighborlist(
            self,
            cutoff=None):
        """ Construct a nearest-neighbor list using ASE:
        https://wiki.fysik.dtu.dk/ase/ase/neighborlist.html

        This function is intended for adaptation with machine learning
        in the future.

        Parameters:
          cutoff: ndarray (96,) or None
            cutoff radius for each of the first 96 elements.

        Returns:
          nl: ASE nearest-neighbor list
            nearest-neighbor lists.
        """

        if cutoff is None:
            # Reduce the radius of sulfur
            r[16] *= 0.80  # S
            r[1] *= 0.80  # H
            cutoff = r[self.atoms.get_atomic_numbers()]

        # Build a nearest neighbor list
        nl = NL(
            cutoff,
            skin=0.2,
            self_interaction=False,
            bothways=False)
        nl.build(self.atoms)

        self.neighbor_list = nl

        return nl

    def id_molecules(
            self,
            cutoff=None):
        """ Identify adsorbed molecules in a given ase atoms object.

        Assumptions:
        - Works via the covalent radius of each atom. When the cutoff overlap,
        a bond is assumed.
        - No atoms with atomic number greater than 21 are adsorbed species.
        - Surface is made of atoms with atomic number > 21 (transition metals).

        WARNING: Current implementation does not support oxides, carbides,
        sulfides, etc...

        Parameters:
          cutoff: ndarray (96,) or None
            cutoff radius for each of the first 96 elements.

        Returns:
          molecules: list (N,)
            List of netowrkx Graph objects representing 2D molecules.
        """

        atoms = self.atoms

        if self.neighbor_list is None:
            self.neighbor_list = self._build_neighborlist(cutoff)

        G = nx.Graph()

        # Find molecules preasent at the surface
        for atom in atoms:

            # Assumes no adsorbates have atomic numbers greater than 21
            if an[atom.symbol] >= 21:
                continue

            a = atom.index
            neighbor_atoms = self.neighbor_list.get_neighbors(a)

            G.add_node(a, symbol=atoms[a].symbol)

            for n in neighbor_atoms[0]:

                # Assumes no adsorbates have atomic numbers greater than 21
                if an[atoms[n].symbol] >= 21:
                    continue

                G.add_edge(a, n)

        molecules = list(nx.connected_component_subgraphs(G))
        self.molecules = molecules

        return molecules

    def id_sites(self):
        """ Estimates the active site of adsorbed molecules.

        Active sites are a reduced order description of an atoms local
        adsorption environment.

        Assumptions:
        - Same as molecule search

        Returns:
          sites: dict
            Keys of molecules found from the molecule search and values
            of dictionaries with index of the molecule atom, index of the
            neighboring metal, and chemical symbol of the neighbor.
        """

        atoms = self.atoms

        molecules = self.molecules
        if self.molecules is None:
            molecules = self.id_molecules()

        # Now attempt to find site descriptions
        sites = {}
        for i, molecule in enumerate(molecules):
            sites[i] = {}

            ids = molecule.nodes()
            neighbors = get_neighbors(atoms, ids)

            for j, nn in neighbors.items():
                sites[i][j] = {}
                for n in nn:
                    if n in ids:
                        continue
                    sites[i][j][n] = atoms[n].symbol

        return sites
