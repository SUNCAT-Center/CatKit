from catkit.flow import qeio
from ase.build import bulk
import unittest
import shutil
import os


class TestQEIO(unittest.TestCase):
    """Test features of catkit.flow.qeio"""

    def tearDown(self):
        """Clear the temp file."""
        if os.path.exists('temp'):
            shutil.rmtree('temp')

    def test_cd(self):
        """Test the cd context manager."""
        with qeio.cd('temp'):
            assert(os.getcwd().split('/')[-1] == 'temp')

    def test_geometry_hash(self):
        """Test the geometry_hash function."""
        atoms = bulk('Pd', cubic=True)
        ghash = qeio.geometry_hash(atoms)

        assert(ghash == '2980adf7416f0035b803eabeef085ecf')


if __name__ == '__main__':
    unittest.main()
