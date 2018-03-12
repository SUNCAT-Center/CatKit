import matplotlib
matplotlib.use('Agg')
from catkit.api.rd_kit import plot_molecule
from catkit.api.rd_kit import get_uff_coordinates
from catkit.pathways import ReactionNetwork
import os


def test_molecule_generation():
    os.makedirs('temp')
    db_name = 'temp/temp.db'
    with ReactionNetwork(db_name=db_name) as rn:
        rn.molecule_search(
            element_pool={'C': 2, 'H': 6},
            multiple_bond_search=False)

        molecules = rn.load_molecules()
        assert (len(molecules) == 17)

        rn.path_search(reconfiguration=False, substitution=False)
        pathways = rn.load_pathways()
        assert (len(pathways) == 27)

        for i, molecule in molecules.items():
            plot_molecule(
                molecule, file_name='temp/molecule-{}.png'.format(i))

            molecule = get_uff_coordinates(molecule, steps=50)
            rn.save_3d_structure(molecule)

        rn.load_3d_structures()
        rn.plot_reaction_network(file_name='temp/reaction-network.png')

    os.unlink(db_name)
    with ReactionNetwork(db_name=db_name) as rn:
        rn.molecule_search(
            element_pool={'C': 2, 'H': 6},
            multiple_bond_search=False)

        molecules = rn.load_molecules()
        assert (len(molecules) == 17)

        rn.path_search(reconfiguration=True, substitution=False)
        pathways = rn.load_pathways()
        assert (len(pathways) == 30)

    os.unlink(db_name)
    with ReactionNetwork(db_name=db_name) as rn:
        rn.molecule_search(
            element_pool={'C': 2, 'H': 6},
            multiple_bond_search=False)

        molecules = rn.load_molecules()
        assert (len(molecules) == 17)

        rn.path_search(reconfiguration=True, substitution=True)
        pathways = rn.load_pathways()
        assert (len(pathways) == 240)

    os.unlink(db_name)
    with ReactionNetwork(db_name=db_name) as rn:
        rn.molecule_search(
            element_pool={'C': 2, 'H': 6},
            multiple_bond_search=True)

        molecules = rn.load_molecules()
        print(len(molecules))
        assert (len(molecules) == 26)

        rn.path_search(reconfiguration=False, substitution=False)
        pathways = rn.load_pathways()
        assert (len(pathways) == 44)

    os.unlink(db_name)
    with ReactionNetwork(db_name=db_name) as rn:
        rn.molecule_search(
            element_pool={'C': 2, 'H': 6},
            multiple_bond_search=True)

        molecules = rn.load_molecules()
        assert (len(molecules) == 26)

        rn.path_search(reconfiguration=True, substitution=False)
        pathways = rn.load_pathways()
        assert (len(pathways) == 48)

    os.unlink(db_name)
    with ReactionNetwork(db_name=db_name) as rn:
        rn.molecule_search(
            element_pool={'C': 2, 'H': 6},
            multiple_bond_search=True)

        molecules = rn.load_molecules()
        assert (len(molecules) == 26)

        rn.path_search(reconfiguration=True, substitution=True)
        pathways = rn.load_pathways()
        assert (len(pathways) == 437)

if __name__ == "__main__":
    test_molecule_generation()
