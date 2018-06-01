from catkit.flow import hpcio
import unittest
import os


class TestHPCIO(unittest.TestCase):
    """Test features of catkit.flow.hpcio"""

    def test_get_server(self):
        """Test get_server function."""
        clusters = ['sherlock', 'slac']

        # Test user defined variables.
        for cluster in clusters:
            os.environ['CLUSTER'] = cluster
            server, local = hpcio.get_server()

            assert(server == cluster)
            del os.environ['CLUSTER']

        globalv = {'LSB_EXEC_CLUSTER': 'slac',
                   'SLURM_CLUSTER_NAME': 'sherlock'}

        # Test user server variables.
        for cluster, name in globalv.items():
            os.environ[cluster] = name
            server, local = hpcio.get_server()

            assert(server == name)
            del os.environ[cluster]

    def assertRaises(ValueError):
        """Test get_server function for failure."""
        hpcio.get_server()

    def test_get_nnodes(self):
        """Test get_nnodes function."""
        servers = ['slac', 'sherlock', 'nersc']

        # Test defines servers
        for server in servers:
            hpcio.get_nnodes(server)

        # Try getting from local environment
        os.environ['CLUSTER'] = 'slac'
        hpcio.get_nnodes()


if __name__ == '__main__':
    unittest.main()
