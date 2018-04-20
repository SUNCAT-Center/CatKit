from catkit import Gratoms
from catgen.utils import to_gratoms
from ase.constraints import FixAtoms
from ase.build import molecule
import networkx as nx


def test_gratoms():
    edges = [(0, 1), (0, 2)]
    atoms = Gratoms(edges=edges)

    for n in atoms.edges:
        assert(n in edges)

    mol = molecule('H2O')
    atoms = to_gratoms(mol)
    atoms.graph.add_edges_from([(0, 1), (0, 2)])

    sym_test = atoms.get_neighbor_symbols(0)
    assert(sym_test.tolist() == ['H', 'H'])

    test_tags = atoms.get_chemical_tags(rank=1)
    assert(test_tags == '2,0,0,0,0,0,0,1')

    test_comp, test_bonds = atoms.get_chemical_tags()
    assert(test_comp == '2,0,0,0,0,0,0,1')
    assert(test_bonds == '4,0,0,0,0,0,0,3')

    atoms.set_constraint(FixAtoms(indices=[0]))
    del atoms[2]
    assert(len(atoms) == 2)

    nx.set_node_attributes(
        atoms.graph, name='valence', values={0: 1, 1: 0})

    test_nodes = atoms.get_unsaturated_nodes(screen=1)
    assert(test_nodes == [0])


if __name__ == "__main__":
    test_gratoms()
