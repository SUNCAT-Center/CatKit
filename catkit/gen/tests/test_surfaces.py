import unittest
from catkit.gen.surface import SlabGenerator
from ase.build import bulk
import numpy as np


class TestSlabGenerator(unittest.TestCase):
    """Test features of catkit.gen.surface"""

    def test_slab_generator(self):
        """Test that no exception is raised when running SlabGenerator."""
        atoms = bulk('Pd', 'fcc', cubic=True)

        SlabGenerator(atoms, miller_index=(1, 1, 1), layers=6, vacuum=4)
        SlabGenerator(atoms, miller_index=(1, 1, 0), layers=6, vacuum=4)
        SlabGenerator(atoms, miller_index=(1, 0, 0), layers=6, vacuum=4)
        SlabGenerator(atoms, miller_index=(0, 0, 0, 1), layers=6)

    def test_terminations(self):
        """Test termination finding on fcc(211)."""
        atoms = bulk('Pd', 'fcc', a=4, cubic=True)
        atoms[3].symbol = 'Cu'

        gen = SlabGenerator(
            atoms, miller_index=(2, 1, 1), layers=8, vacuum=4)

        gen.get_slab(iterm=1)

    def test_gcd_miller_index(self):
        """Test that miller indices are properly reduces to their
        greatest common divisor.
        """
        atoms = bulk('Pd', 'fcc')
        gen = SlabGenerator(atoms, miller_index=(3, 3, 3),
                            layers=3, layer_type='trim')

        np.testing.assert_array_equal(gen.miller_index, [1, 1, 1])

    def test_no_standardization(self):
        """Test that the slab generator functions if not standardized."""
        atoms = bulk('Pd', 'fcc', cubic=True)
        SlabGenerator(atoms, miller_index=(1, 1, 1), layer_type='trim',
                      layers=3, standardize_bulk=True)

    def test_matrix_notation_search(self):
        """Test the matrix notation algorithm."""
        atoms = bulk('Pd', 'fcc', cubic=True)
        gen = SlabGenerator(atoms, miller_index=(1, 1, 1), layers=6)

        gen.get_slab(size=1)
        gen.get_slab(size=[[1, 0], [0, 1]])


if __name__ == '__main__':
    unittest.main()
