import unittest
from catkit.build import molecule
from catgen.build import surface
from ase.build import bulk
import numpy as np
from np.testing import assert_array_equal


class TestBuild(unittest.TestCase):
    """Test functions in the build module."""

    def test_surface(self):
        """Test the helper surface generator."""
        slab = surface('Pd', size=(2, 2, 4), vacuum=10)

        # Slab should have 16 Pd atoms
        assert(len(slab) == 16)

        correct_connectivity = np.array([
            [0, 2, 2, 2, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [2, 0, 2, 2, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [2, 2, 0, 2, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [2, 2, 2, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 1, 1, 0, 2, 2, 2, 1, 1, 1, 0, 0, 0, 0, 0],
            [1, 0, 1, 1, 2, 0, 2, 2, 1, 1, 0, 1, 0, 0, 0, 0],
            [1, 1, 0, 1, 2, 2, 0, 2, 1, 0, 1, 1, 0, 0, 0, 0],
            [1, 1, 1, 0, 2, 2, 2, 0, 0, 1, 1, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 1, 1, 0, 0, 2, 2, 2, 1, 1, 1, 0],
            [0, 0, 0, 0, 1, 1, 0, 1, 2, 0, 2, 2, 1, 1, 0, 1],
            [0, 0, 0, 0, 1, 0, 1, 1, 2, 2, 0, 2, 1, 0, 1, 1],
            [0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 2, 2, 2],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 2, 0, 2, 2],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 2, 2, 0, 2],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 0]])
        assert_array_equal(slab.connectivity, correct_connectivity)

        correct_surf_atoms = np.array([12, 13, 14, 15])
        assert_array_equal(slab.get_surface_atoms(), correct_surf_atoms)

        # Test the ability to pass an atoms object
        atoms = bulk('Pd', cubic=True)
        slab = surface(atoms, size=(2, 2, 4), vacuum=10)

    def test_molecule(self):
        """Test the helper molecule generator."""
        images = molecule('H')
        assert(len(images) == 1)

        images = molecule('C3NH16', vacuum=10)
        assert(len(images) == 4)


if __name__ == '__main__':
    unittest.main()
