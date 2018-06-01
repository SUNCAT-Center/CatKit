import unittest
from catkit.flow.qeio import cd
import shutil
import os


class TestQEIO(unittest.TestCase):
    """Test features of catkit.flow.qeio"""

    def tearDown(self):
        """Clear the temp file."""
        shutil.rmtree('temp')

    def test_cd(self):
        """Test the cd context manager."""
        with cd('temp'):
            assert(os.getcwd().split('/')[-1] == 'temp')


if __name__ == '__main__':
    unittest.main()
