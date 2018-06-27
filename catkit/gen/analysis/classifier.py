from catkit import Gratoms
from .. import utils
from ase import Atoms
import networkx as nx
import numpy as np


class Classifier():
    """Class for classification of various aspects of an an atomic
    unit cell.

    Currently, a tool for classification of adsorbates on surface
    environments and the active sites they rest on.
    """

    def __init__(self, atoms):
        """Return unique coordinate values of a given atoms object
        for a specified axis.

        Parameters
        ----------
            atoms : atoms object
        """
        self.atoms = atoms
        self.ads_atoms = None
        self.slab_atoms = None
        self.surface_atoms = None

    def id_slab_atoms(
            self,
            classifier='trivial',
            tag=False,
            rtol=1e-3):
        """Return the indices of the slab atoms using select characterization
        techniques.

        Parameters
        ----------
        classifier : str
            Classification technique to identify slab atoms.

            'trivial':
            Slab atoms assumed to have atomic number == 13 or >= 21.
        tag : bool
            Return adsorbate atoms with tags of 2.
        rtol : float
            Relative cutoff distance for tagging layers.

        Returns
        -------
        slab_atoms : ndarray (n,)
            Index of slab atoms found.
        """
        atoms = self.atoms

        if classifier == 'trivial':
            slab_atoms = np.where((atoms.numbers == 13) |
                                  (atoms.numbers >= 21))[0]

        if tag:
            zpos = np.sort(atoms.positions[slab_atoms][:, -1])
            new_tags = np.zeros_like(zpos, dtype=int)
            tag = 1
            for i, z in enumerate(zpos):
                if new_tags[i] != 0:
                    continue
                layer = np.isclose(z, zpos, rtol=rtol)
                new_tags[layer] = tag
                tag += 1

            tags = self.atoms.get_tags()
            tags[slab_atoms] = new_tags[::-1]
            self.atoms.set_tags(tags)

        self.slab_atoms = slab_atoms

        return slab_atoms

    def id_adsorbate_atoms(self, classifier='trivial', tag=False):
        """Identify adsorbed atoms in a given atoms object.

        Parameters
        ----------
        classifier : str
            Classification technique to identify adsorbate atoms.

            'trivial':
            Adsorbate atoms assumed to have atomic number != 13 or < 21.
        tag : bool
           Return adsorbate atoms with tags of -2.

        Returns
        -------
        ads_atoms : ndarray (n,)
            Index of adsorbate atoms found.
        """
        atoms = self.atoms

        if classifier == 'trivial':
            ads_atoms = np.where((atoms.numbers != 13) &
                                 (atoms.numbers < 21))[0]

        if tag:
            tags = self.atoms.get_tags()
            tags[ads_atoms] = -2
            self.atoms.set_tags(tags)

        self.ads_atoms = ads_atoms

        return ads_atoms

    def id_surface_atoms(self, classifier='voronoi_sweep'):
        """Identify surface atoms of an atoms object. This will
        require that adsorbate atoms have already been identified.

        Parameters
        ----------
        classifier : str
            Classification technique to identify surface atoms.

            'voronoi_sweep':
            Create a sweep of proxy atoms above surface. Surface atoms
            are those which are most frequent neighbors of the sweep.

        Returns
        -------
        surface_atoms : ndarray (n,)
            Index of the surface atoms in the object.
        """
        atoms = self.atoms.copy()

        # Remove adsorbates before analysis
        ads_atoms = self.ads_atoms
        if ads_atoms is None:
            ads_atoms = self.id_adsorbate_atoms()
        del atoms[ads_atoms]

        if classifier == 'voronoi_sweep':
            spos = atoms.get_scaled_positions()
            zmax = np.max(spos[:, -1])

            # Create a distribution of points to screen with
            # 2.5 angstrom defines the absolute separation
            dvec = (np.linalg.norm(atoms.cell[:-1], axis=1) / 2.5) ** -1
            xy = np.mgrid[0:1:dvec[0], 0:1:dvec[1]].reshape(2, -1)
            z = np.ones_like(xy[0]) * zmax
            xyz = np.vstack((xy, z)).T

            screen = np.dot(xyz, atoms.cell)

            n = len(atoms)
            m = len(screen)
            ind = np.arange(n, n + m)

            slab_atoms = np.arange(n)

            satoms = []
            # 2 - 3 Angstroms seems to work for a large range of indices.
            for k in np.linspace(2, 3, 10):
                wall = screen.copy() + [0, 0, k]

                atm = Atoms(['X'] * m, positions=wall)
                test_atoms = atoms + atm

                con = utils.get_voronoi_neighbors(test_atoms)
                surf_atoms = np.where(con[ind].sum(axis=0)[slab_atoms])[0]
                satoms += [surf_atoms]

            len_surf_atoms = [len(_) for _ in satoms]
            uni, ind, cnt = np.unique(
                len_surf_atoms, return_counts=True, return_index=True)

            max_cnt = np.argmax(cnt)
            surf_atoms = satoms[ind[max_cnt]]

        self.surface_atoms = surf_atoms

        return surf_atoms

    def id_adsorbates(self, classifier='radial', return_atoms=False):
        """Return a list of Gratoms objects for each adsorbate
        classified on a surface. Requires classification of adsorbate
        atoms.

        Parameters
        ----------
        classifier : str
            Classification technique to identify individual adsorbates.

            'radial':
            Use standard cutoff distances to identify neighboring atoms.

        return_atoms : bool
            Return Gratoms objects instead of adsorbate indices.

        Returns
        -------
        adsorabtes : list (n,)
            Adsorbate indices of adsorbates in unit cell.
        """
        atoms = self.atoms.copy()

        # Remove the slab atoms
        ads_atoms = self.ads_atoms
        if ads_atoms is None:
            ads_atoms = self.id_adsorbate_atoms()

        if classifier == 'radial':
            con = utils.get_cutoff_neighbors(atoms)
            ads_con = con[ads_atoms][:, ads_atoms]
            G = nx.Graph()
            G.add_nodes_from(ads_atoms)
            edges = utils.connectivity_to_edges(ads_con, indices=ads_atoms)
            G.add_weighted_edges_from(edges, weight='bonds')
            SG = nx.connected_component_subgraphs(G)

            adsorabtes = []
            for sg in SG:
                nodes = list(sg.nodes)
                if return_atoms:
                    edges = list(sg.edges)
                    ads = Gratoms(
                        numbers=atoms.numbers[nodes],
                        positions=atoms.positions[nodes],
                        edges=edges)
                    ads.center(vacuum=5)
                else:
                    ads = nodes
                adsorabtes += [ads]

        return adsorabtes
