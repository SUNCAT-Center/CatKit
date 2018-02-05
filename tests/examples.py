from catkit.surface import SlabGenerator
from catkit.pathways import ReactionNetwork
from ase.build import bulk
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
    )

    expected_number = {'top': 2, 'bridge': 3, 'hollow': 4, '4fold': 0}

    for k, v in sites.items():
        positions, points, _ = v
        assert(len(positions) == expected_number[k])


def test_molecule_generation():
    db_name = 'temp.db'

    with ReactionNetwork(db_name=db_name) as rn:

        molecules = rn.molecule_search(
            element_pool={'C': 2, 'H': 6},
            multiple_bond_search=False)
        rn.save_molecules(molecules)

        molecules = rn.load_molecules()

        assert(len(molecules) == 17)

        pathways = rn.path_search(
            reconfiguration=True,
            substitution=True)

        pathways = rn.load_pathways()

        assert(len(pathways) == 240)

    os.unlink('temp.db')


if __name__ == "__main__":
    test_slab_generation()
    test_adsorption_sites()
    test_molecule_generation()
