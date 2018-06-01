import unittest
from catkit.gen.pathways import ReactionNetwork
import shutil
import os

db_name = 'temp/temp.db'


class TestSlabGenerator(unittest.TestCase):
    """Test features of catkit.gen.pathways"""

    def setUp(self):
        """Set up a temp file."""
        if not os.path.exists('temp'):
            os.makedirs('temp')

    def tearDown(self):
        """Clear the temp file."""
        shutil.rmtree('temp')

    def test_isomorphic_tree_generation(self):
        """Test simple isomorphic-tree search."""
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

            rn.plot_reaction_network(file_name='temp/reaction-network.png')
        os.unlink(db_name)

    def test_advanced_paths(self):
        """Test reconfiguration and substitution step enumeration."""
        with ReactionNetwork(db_name=db_name) as rn:
            rn.molecule_search(
                element_pool={
                    'C': 2,
                    'H': 6
                }, multiple_bond_search=False)

            molecules = rn.load_molecules()
            assert (len(molecules) == 17)

            rn.path_search(reconfiguration=True, substitution=True)
            pathways = rn.load_pathways()
            assert (len(pathways) == 240)
        os.unlink(db_name)

    def test_advanced_paths_wloops(self):
        """Test reconfiguration and substitution step enumeration with
        multiple bonding molecules.
        """
        with ReactionNetwork(db_name=db_name) as rn:
            rn.molecule_search(
                element_pool={
                    'C': 2,
                    'H': 6
                }, multiple_bond_search=True)

            molecules = rn.load_molecules()
            assert (len(molecules) == 26)

            rn.path_search(reconfiguration=True, substitution=True)
            pathways = rn.load_pathways()
            assert (len(pathways) == 437)


if __name__ == '__main__':
    unittest.main()
