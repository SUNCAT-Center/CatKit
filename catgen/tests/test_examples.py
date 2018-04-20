import matplotlib
matplotlib.use('Agg')
from catgen.surface import SlabGenerator
from catgen.adsorption import AdsorptionSites
from catgen.pathways import ReactionNetwork
from catgen.api.rd_kit import plot_molecule
from catgen.api.rd_kit import get_uff_coordinates
from ase.utils import formula_hill
from ase.build import bulk
import networkx as nx
import numpy as np
import shutil
import os


def test_surface_examples():
    atoms = bulk('Pd', 'fcc', a=4, cubic=True)
    atoms[3].symbol = 'Cu'

    gen = SlabGenerator(
        atoms, miller_index=(2, 1, 1), layers=9, fixed=5, vacuum=4)

    terminations = gen.get_unique_terminations()
    assert (len(terminations) == 2)

    for i, t in enumerate(terminations):
        slab = gen.get_slab(iterm=i)
        assert (len(slab) == 18)

    atoms = bulk('Pd', 'hcp', a=3, cubic=True)

    gen = SlabGenerator(
        atoms, miller_index=(1, 1, 0), layers=6, fixed=2, vacuum=4)

    atoms = gen.get_slab()
    con_matrix = gen.get_graph_from_bulk(atoms, attach=True)

    test_con_matrix = np.array(
        [[0.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0], [2.0, 0.0, 2.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
         [2.0, 2.0, 0.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0], [2.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
         [1.0, 0.0, 2.0, 2.0, 0.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0,
          0.0], [0.0, 1.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 1.0, 0.0, 0.0], [
              0.0, 0.0, 1.0, 0.0, 2.0, 2.0, 0.0, 2.0, 2.0, 2.0, 1.0, 0.0
          ], [0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0, 0.0, 1.0], [
              0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 2.0, 2.0, 0.0, 2.0, 2.0, 2.0
          ], [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0], [
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 2.0, 2.0, 0.0, 2.0
          ], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 0.0]])

    np.testing.assert_allclose(con_matrix, test_con_matrix)

    # We can identify both top and bottom sites.
    top, bottom = gen.get_voronoi_surface_atoms(atoms)
    atoms.set_surface_atoms(top)

    test_surf_atoms = np.array([8, 9, 10, 11])
    np.testing.assert_allclose(top, test_surf_atoms)

    atoms = bulk('Pd', 'fcc', a=5, cubic=True)
    atoms[3].symbol = 'Cu'

    gen = SlabGenerator(
        atoms, miller_index=(1, 1, 1), layers=3, fixed=2, vacuum=10)

    atoms = gen.get_slab(primitive=True)
    coordinates, connectivity = gen.adsorption_sites(atoms)

    test_connectivity = np.array([1, 1, 2, 2, 2, 3, 3, 3, 3])
    np.testing.assert_allclose(connectivity, test_connectivity)


def test_adsorption_examples():
    atoms = bulk('Pd', 'fcc', a=5, cubic=True)
    atoms[3].symbol = 'Cu'

    gen = SlabGenerator(atoms, miller_index=(1, 1, 1), layers=3, vacuum=4)

    atoms = gen.get_slab(primitive=True)
    atoms.set_surface_atoms([8, 9, 10, 11])

    sites = AdsorptionSites(atoms)
    sites.plot('./images/Pd3Cu-adsorption-sites.png')

    atoms = bulk('Pd', 'fcc', a=5, cubic=True)
    atoms[3].symbol = 'Cu'

    gen = SlabGenerator(atoms, miller_index=(3, 2, 1), layers=13, vacuum=5)

    atoms = gen.get_slab(primitive=True)

    top, _ = gen.get_voronoi_surface_atoms(atoms)
    atoms.set_surface_atoms(top)
    sites = AdsorptionSites(atoms)

    coordinates = sites.get_coordinates()
    assert (len(coordinates) == 56)
    connectivity = sites.get_connectivity()
    test_connectivity = np.array([
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 4, 4
    ])
    np.testing.assert_allclose(connectivity, test_connectivity)

    topology = sites.get_topology()
    assert (len(topology) == 56)

    periodic = sites.get_periodic_sites()
    symmetric = sites.get_symmetric_sites()
    print(periodic)
    np.testing.assert_allclose(symmetric, periodic)

    atoms = bulk('Pd', 'fcc', a=5, cubic=True)
    atoms[3].symbol = 'Cu'

    gen = SlabGenerator(atoms, miller_index=(2, 1, 1), layers=10, vacuum=5)

    atoms = gen.get_slab(primitive=True)

    top, _ = gen.get_voronoi_surface_atoms(atoms, attach_graph=False)
    atoms.set_surface_atoms(top)
    sites = AdsorptionSites(atoms)

    coordinates = sites.get_coordinates()
    vectors = sites.get_adsorption_vectors(unique=False)
    assert (len(vectors) == 32)


def test_gas_phase_example():
    db_name = 'C2H6-example.db'
    with ReactionNetwork(db_name=db_name) as rn:

        rn.molecule_search(
            element_pool={
                'C': 2,
                'H': 6
            }, multiple_bond_search=False)
        molecules = rn.load_molecules()

        assert (len(molecules) == 17)

        for i, molecule in molecules.items():
            plot_molecule(
                molecule, file_name='./images/molecule-{}.png'.format(i))

            molecule = get_uff_coordinates(molecule, steps=50)
            rn.save_3d_structure(molecule)

        images = rn.load_3d_structures()

        assert (len(images) == 17)

        rn.path_search(reconfiguration=False, substitution=False)
        rn.plot_reaction_network(file_name='./images/reaction-network.png')
        pathways = rn.load_pathways()

        assert (len(pathways) == 27)

    names = np.empty(len(molecules) + 1, dtype='a5')
    names[0] = ''
    for k, v in molecules.items():
        atn = nx.get_node_attributes(v.graph, 'number')

        hill = formula_hill(list(atn.values()))
        names[k] = hill

    for path in pathways:
        print('|{} + {} --> {} + {}|'.format(*names[path]))

    os.unlink(db_name)


if __name__ == "__main__":
    test_surface_examples()
    test_adsorption_examples()
    test_gas_phase_example()
