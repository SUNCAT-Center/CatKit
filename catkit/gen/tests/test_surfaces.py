import unittest
from catkit.gen.surface import SlabGenerator
from ase.build import bulk


class TestSlabGenerator(unittest.TestCase):
    """Test features of catkit.gen.surface"""

    def test_slab_generator(self):
        """Test that no exception is raised when running SlabGenerator."""
        atoms = bulk('Pd', 'fcc', a=4, cubic=True)

        SlabGenerator(atoms, miller_index=(1, 1, 1), layers=6, vacuum=4)
        SlabGenerator(atoms, miller_index=(1, 1, 0), layers=6, vacuum=4)
        SlabGenerator(atoms, miller_index=(1, 0, 0), layers=6, vacuum=4)

    def test_terminations(self):
        """Test termination finding on fcc(211)."""
        atoms = bulk('Pd', 'fcc', a=4, cubic=True)
        atoms[3].symbol = 'Cu'

        gen = SlabGenerator(
            atoms, miller_index=(2, 1, 1), layers=8, vacuum=4)

        gen.get_slab(iterm=1)


if __name__ == '__main__':
    unittest.main()
