from catkit.db import FingerprintDB
from ase.db import connect
import networkx as nx
import numpy as np
import json

db = connect('data/molecules.db')
with open('data/chemical-properties.json', 'r') as f:
    properties = json.load(f)


def test_database():
    dstar = 8
    with FingerprintDB('fingerprints.db') as fpd:
        for i in range(dstar + 1):
            fpd.parameter_entry('Z{}'.format(i),
                                'Atomic charge: depth {}'.format(i))
            fpd.parameter_entry('r{}'.format(i),
                                'Cordero covalent radius: depth {}'.format(i))
            fpd.parameter_entry('x{}'.format(i),
                                'Pauling electronegetivity: depth {}'.format(i))
            fpd.parameter_entry('T{}'.format(i),
                                'Coordination number: depth {}'.format(i))
            fpd.parameter_entry('1{}'.format(i),
                                'Unity: depth {}'.format(i))
        fpd.parameter_entry('Ef', 'Formation energy')
        fpd.parameter_entry('Et', 'Total energy')

        par = fpd.get_parameters()

        for d in db.select('id<100'):
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
                ind = str(int(n))
                W[1][i] = properties[ind]['covalent_radius_cordero'] / 100
                W[2][i] = properties[ind]['en_pauling']
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


if __name__ == "__main__":
    test_database()
