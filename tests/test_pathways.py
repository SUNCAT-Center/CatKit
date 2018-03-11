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
            element_pool={
                'C': 2,
                'H': 6
            }, multiple_bond_search=True)

        molecules = rn.load_molecules()
        assert (len(molecules) == 26)

        rn.path_search(reconfiguration=True, substitution=True)

        pathways = rn.load_pathways()
        # The appropriate length of the this needs verification.
        assert (len(pathways) == 618)

        for i, molecule in molecules.items():
            plot_molecule(
                molecule, file_name='temp/molecule-{}.png'.format(i))

            molecule = get_uff_coordinates(molecule, steps=50)
            rn.save_3d_structure(molecule)

        rn.load_3d_structures()
        rn.plot_reaction_network(file_name='temp/reaction-network.png')


if __name__ == "__main__":
    test_molecule_generation()
