import unittest
from catkit.gen.analysis import Classifier
from catkit.build import molecule


class TestAnalysis(unittest.TestCase):
    """Test features of the analysis module."""

    def test_classifier():
        """Test identification of adsorbate atoms and slab atoms."""
        atoms = molecule('C2H2')[0]

        # This needs to be made rigorous.
        cl = Classifier(atoms)
        cl.id_adsorbate_atoms()
        cl.id_slab_atoms()


if __name__ == '__main__':
    unittest.main()
