import unittest
from catkit.build import molecule
from catkit.build import surface
from ase.build import bulk
import numpy as np
from numpy.testing import assert_array_equal


class TestBuild(unittest.TestCase):
    """Test functions in the build module."""

    def test_surface(self):
        """Test catkit.build.surface generator."""
        slab = surface('Pd', size=(2, 2, 4), vacuum=10)

        # Slab should have 16 Pd atoms
        assert(len(slab) == 16)
        degree = slab.degree

        test_degree = np.array(
            [9, 9, 9, 9, 12, 12, 12, 12, 12, 12, 12, 12, 9, 9, 9, 9])

        assert_array_equal(degree, test_degree)

        correct_surf_atoms = np.array([12, 13, 14, 15])
        assert_array_equal(slab.get_surface_atoms(), correct_surf_atoms)

        # Test the ability to pass an atoms object
        atoms = bulk('Pd', cubic=True)
        slab = surface(atoms, size=(2, 2, 4), vacuum=10)

        # Test orthogonalization search
        slab = surface(atoms, size=(1, 4), vacuum=10)

    def test_molecule(self):
        """Test catkit.build.molecule generator."""
        images = molecule('H')
        assert(len(images) == 1)

        images = molecule('C3NH16', vacuum=10)
        assert(len(images) == 4)


if __name__ == '__main__':
    unittest.main()
