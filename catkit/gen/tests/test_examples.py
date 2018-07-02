import unittest
from catkit.build import surface
from catkit.gen.surface import SlabGenerator
from catkit.gen.adsorption import AdsorptionSites
from catkit.gen.pathways import ReactionNetwork
from ase.utils import formula_hill
from ase.build import bulk
import networkx as nx
import numpy as np
from glob import glob
import os


class TestGenDocs(unittest.TestCase):
    """Test the examples in the documentation."""

    def tearDown(self):
        """Remove any lingering files"""
        for f in glob('*.db'):
            os.unlink(f)

        for f in glob('*.png'):
            os.unlink(f)

    def test_surface_examples(self):
        """Test surface construction examples."""
        atoms = bulk('Pd', 'fcc', a=4, cubic=True)
        atoms[3].symbol = 'Cu'

        gen = SlabGenerator(
            atoms, miller_index=(2, 1, 1), layers=6, fixed=5, vacuum=4)

        terminations = gen.get_unique_terminations()
        assert (len(terminations) == 2)

        for i, t in enumerate(terminations):
            slab = gen.get_slab(iterm=i)
            assert (len(slab) == 16)

        atoms = surface('Pd', size=(3, 3), miller=(1, 1, 1), vacuum=4)
        con_matrix = atoms.connectivity

        test_con_matrix = np.array([
            [0, 3, 3, 1, 1, 1, 0, 0, 0],
            [3, 0, 3, 1, 1, 1, 0, 0, 0],
            [3, 3, 0, 1, 1, 1, 0, 0, 0],
            [1, 1, 1, 0, 3, 3, 1, 1, 1],
            [1, 1, 1, 3, 0, 3, 1, 1, 1],
            [1, 1, 1, 3, 3, 0, 1, 1, 1],
            [0, 0, 0, 1, 1, 1, 0, 3, 3],
            [0, 0, 0, 1, 1, 1, 3, 0, 3],
            [0, 0, 0, 1, 1, 1, 3, 3, 0]])

        np.testing.assert_allclose(con_matrix, test_con_matrix)

        test_surf_atoms = np.array([6, 7, 8])
        np.testing.assert_allclose(atoms.get_surface_atoms(), test_surf_atoms)

        atoms = bulk('Pd', 'fcc', a=5, cubic=True)
        atoms[3].symbol = 'Cu'

        gen = SlabGenerator(
            atoms, miller_index=(1, 1, 1), layers=3, layer_type='trim',
            fixed=2, vacuum=10)

        atoms = gen.get_slab()
        coordinates, connectivity = gen.adsorption_sites(atoms)

        test_connectivity = np.array([1, 1, 2, 2, 2, 3, 3, 3, 3])
        np.testing.assert_allclose(connectivity, test_connectivity)

    def test_adsorption_examples(self):
        """Test adsorption structure examples."""
        atoms = bulk('Pd', 'fcc', a=5, cubic=True)
        atoms[3].symbol = 'Cu'

        gen = SlabGenerator(atoms, miller_index=(1, 1, 1), layer_type='trim',
                            layers=3, vacuum=4)

        atoms = gen.get_slab()

        sites = AdsorptionSites(atoms)
        sites.plot('./Pd3Cu-adsorption-sites.png')

        atoms = bulk('Pd', 'fcc', a=5, cubic=True)
        atoms[3].symbol = 'Cu'

        gen = SlabGenerator(atoms, miller_index=(3, 2, 1), layers=8, vacuum=5)

        atoms = gen.get_slab()
        sites = AdsorptionSites(atoms)

        coordinates = sites.get_coordinates()
        assert (len(coordinates) == 56)
        connectivity = sites.get_connectivity()
        test_connectivity = np.array([
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4])
        np.testing.assert_allclose(connectivity, test_connectivity)

        topology = sites.get_topology()
        assert (len(topology) == 56)

        periodic = sites.get_periodic_sites()
        symmetric = sites.get_symmetric_sites()
        np.testing.assert_allclose(symmetric, periodic)

        atoms = bulk('Pd', 'fcc', a=5, cubic=True)
        atoms[3].symbol = 'Cu'

        gen = SlabGenerator(atoms, miller_index=(2, 1, 1), layers=8, vacuum=5)

        atoms = gen.get_slab()
        sites = AdsorptionSites(atoms)

        coordinates = sites.get_coordinates()
        vectors = sites.get_adsorption_vectors(unique=False)
        assert (len(vectors) == 32)

    def test_gas_phase_example(self):
        """Test the gas phase enumerator examples."""
        db_name = 'C2H6-example.db'
        with ReactionNetwork(db_name=db_name) as rn:

            rn.molecule_search(
                element_pool={
                    'C': 2,
                    'H': 6
                }, multiple_bond_search=False)
            molecules = rn.load_molecules()

            assert (len(molecules) == 17)

            rn.path_search(reconfiguration=False, substitution=False)
            pathways = rn.load_pathways()

            assert (len(pathways) == 27)

        names = np.empty(len(molecules) + 1, dtype='a5')
        names[0] = ''
        for k, v in molecules.items():
            atn = nx.get_node_attributes(v.graph, 'number')

            hill = formula_hill(list(atn.values()))
            names[k] = hill

        for path in pathways:
            print('|{} + {} --> {} + {}|'.format(*names[path]))


if __name__ == '__main__':
    unittest.main()
