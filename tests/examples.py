import matplotlib
matplotlib.use('Agg')
from catkit.surface import SlabGenerator
from catkit.pathways import plot_molecule
from catkit.pathways import ReactionNetwork
from ase.build import bulk
import shutil
import os


def test_slab_generation():
    atoms = bulk('Pd', 'fcc', a=4, cubic=True)
    atoms[3].symbol = 'Cu'

    Gen = SlabGenerator(
        atoms,
        miller_index=[2, 1, 1],
        layers=9,
        fixed=5,
        vacuum=10,
    )

    terminations = Gen.get_unique_terminations()

    images = []
    for i, t in enumerate(terminations):
        slab = Gen.get_slab(iterm=i)
        slab.center(axis=2, vacuum=5)

        images += [slab]

    assert(len(images) == 2)


def test_adsorption_sites():
    atoms = bulk('Pd', 'fcc', a=5, cubic=True)
    atoms[3].symbol = 'Cu'

    gen = SlabGenerator(
        atoms,
        miller_index=(1, 1, 1),
        layers=3,
        fixed=2,
        vacuum=10,
    )

    slab = gen.get_slab(
        primitive=True,
    )

    sites = gen.adsorption_sites(
        slab,
        symmetry_reduced=True,
        vectors=True,
    )

    expected_number = {'top': 2, 'bridge': 3, 'hollow': 4, '4fold': 0}

    for k, v in sites.items():
        positions, points, _ = v
        assert(len(positions) == expected_number[k])


def test_molecule_generation():
    os.makedirs('temp')
    db_name = 'temp\temp.db'

    with ReactionNetwork(db_name=db_name) as rn:

        molecules = rn.molecule_search(
            element_pool={'C': 2, 'H': 6},
            multiple_bond_search=True)
        rn.save_molecules(molecules)

        molecules = rn.load_molecules()

        assert(len(molecules) == 26)

        pathways = rn.path_search(
            reconfiguration=True,
            substitution=True)

        pathways = rn.load_pathways()

        assert(len(pathways) == 437)

        molecules = rn.load_molecules()
        for i, molecule in molecules.items():
            rn.save_3d_structure(molecule, uff=50)
            plot_molecule(
                molecule,
                file_name='./temp/molecule-{}.png'.format(i))

            rn.load_3d_structures()

        rn.plot_reaction_network(file_name='./temp/reaction-network.png')

    # Cleanup
    shutil.rmtree('temp')


if __name__ == "__main__":
    test_slab_generation()
    test_adsorption_sites()
    test_molecule_generation()
