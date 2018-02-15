from .utils import get_neighbors
import networkx as nx
import networkx.algorithms.isomorphism as iso
from ase.neighborlist import NeighborList as NL
from ase.data import atomic_numbers as an
from ase.data import covalent_radii as r
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
            if an[atom.symbol] >= 21 or an[atom.symbol] == 13:
                continue

            a = atom.index
            neighbor_atoms = self.neighbor_list.get_neighbors(a)

            G.add_node(
                a,
                number=atom.number,
                symbol=atom.symbol)

            for n in neighbor_atoms[0]:

                # Assumes no adsorbates have atomic numbers greater than 21
                if an[atoms[n].symbol] >= 21 or an[atom.symbol] == 13:
                    continue

                G.add_edge(a, n, bonds=1)

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


def id_reconstruction(
        images,
        rmean=4,
        save=False):
    """ Identify a reconstruction even analyzing changes in the forces.

    Parameters:
      images: list of ASE atoms objects
        Relaxation trajectory.

      rmean: int
        Number of values to use for rolling mean.

      show: bool
        Create a figure to display the events located.

    Returns:
      predicted_events: list of int
        index of images predicted before the event occurs.
    """

    forces = []
    for i, atoms in enumerate(images):
        forces += [np.sqrt((atoms.get_forces() ** 2).sum())]
    forces = np.array(forces)

    frm = pd.rolling_mean(forces, 4)
    fdiff = np.diff(frm)
    fterm = np.array([fdiff > 0.25 * frm[:-1]]).astype(int)[0]
    predicted_events = np.where(fterm[:-1] < fterm[1:])[0]

    if save:
        fig, ax = plt.subplots(figsize=(6, 4))
        l, = plt.plot(range(1, len(images) + 1), frm)
        ax.fill_between(
            range(1, len(images) + 1),
            np.zeros(len(images)),
            frm,
            facecolor=l.get_color(),
            alpha=0.5,
            interpolate=True)
        for i in predicted_events:
            plt.text(i - 1, 0.9, i)
            plt.axvline(i, ls='--', color='0.4')

        ylim = ax.get_ylim()
        plt.xlim(4, len(images))
        plt.ylim(0, ylim[1])
        plt.xlabel('Relaxation step')
        plt.ylabel('Force running mean (eV/$\AA$)')
        plt.savefig(save)
        plt.close()

    return predicted_events


def reactant_indices(R1, R2, P, broken_bond):
    """ Match the indices of a pair of reactants from a
     product xand broken bond.

    Parameters:
      R1: networkx MultiGraph
        Graph representing reactant 1
      R2: networkx MultiGraph
        Graph representing reactant 2
      P: networkx MultiGraph
        Graph representing the product
      broken_bond: list (2,)
        Indices representing the edge of the product
        to be removed.

    Returns: ndarrays (n,)
      Indices of the product graph sorted by the order of
      the reactants indices.
    """

    GM = nx.algorithms.isomorphism.GraphMatcher
    em = iso.numerical_edge_match('bonds', 1)
    nm = iso.numerical_node_match('number', 1)

    Pgraph = P.copy()
    u, v = broken_bond
    Pgraph.graph.remove_edge(u, v)
    Rgraph = R1 + R2

    gm = GM(
        Pgraph.graph,
        Rgraph.graph,
        edge_match=em,
        node_match=nm
    )

    gm.is_isomorphic()

    pindex = np.empty(len(Pgraph), dtype=int)
    for k, v in gm.mapping.items():
        pindex[k] = v

    return pindex
