import unittest
from catkit.pawprint.db import FingerprintDB
from ase.db import connect
import networkx as nx
import numpy as np
import pkg_resources
import json
import os

# Load in properties for testing.
db_path = pkg_resources.resource_filename(
    'catkit', 'tests/data/molecules.db')
db = connect(db_path)

data_path = pkg_resources.resource_filename(
    'catkit', 'data/properties.json')
with open(data_path, 'r') as f:
    properties = json.load(f)


class TestDB(unittest.TestCase):
    """Test features of the learn module."""

    def tearDown(self):
        """Remove the temporary database."""
        os.unlink('tmp-fingerprints.db')

    def test_database(self):
        """Test a general use case for the fingerprint database."""
        dstar = 8
        with FingerprintDB('tmp-fingerprints.db') as fpd:
            for i in range(dstar + 1):
                fpd.parameter_entry('Z{}'.format(i),
                                    'Atomic charge: depth {}'.format(i))
                fpd.parameter_entry(
                    'r{}'.format(i),
                    'Cordero covalent radius: depth {}'.format(i))
                fpd.parameter_entry(
                    'x{}'.format(i),
                    'Pauling electronegetivity: depth {}'.format(i))
                fpd.parameter_entry('T{}'.format(i),
                                    'Coordination number: depth {}'.format(i))
                fpd.parameter_entry(
                    '1{}'.format(i), 'Unity: depth {}'.format(i))
            fpd.parameter_entry('Ef', 'Formation energy')
            fpd.parameter_entry('Et', 'Total energy')

            par = fpd.get_parameters()

            for d in db.select():
                fpd.image_entry(d)

                atoms = d.toatoms()

                edges = [tuple([u, v]) for u, v in d.data['edges']]
                G = nx.Graph()
                G.add_nodes_from(range(len(atoms)))
                G.add_edges_from(edges)
                distance_matrix = nx.floyd_warshall_numpy(G)
                Bm = np.zeros(distance_matrix.shape)

                W = np.ones((5, len(atoms)))
                W[0] = atoms.numbers

                for i, n in enumerate(W[0]):
                    W[1][i] = properties['covalent_radius_cordero'][int(n)]
                    W[2][i] = properties['en_pauling'][int(n)]
                    W[3][i] = len(G[i])

                for dd in range(dstar + 1):
                    B = Bm.copy()
                    B[distance_matrix == dd] = 1
                    AC = np.dot(np.dot(W, B), W.T).diagonal()

                    for j, v in enumerate(AC):
                        ind = j + dd * len(AC)
                        fpd.fingerprint_entry(d.id, par[ind], v)

                fpd.fingerprint_entry(d.id, 46, d.Uref)
                fpd.fingerprint_entry(d.id, 47, atoms.get_potential_energy())


if __name__ == '__main__':
    unittest.main()
