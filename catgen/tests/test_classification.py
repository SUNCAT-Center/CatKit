import matplotlib
matplotlib.use('Agg')
from ase.io import read
from catgen.classification import Classifier


def test_classifier():
    atoms = read('data/molecules.db')

    # This needs to be made rigorous, but also to be converted to
    # Gratoms automatically.
    ident = Classifier(atoms)
    mol = ident.id_molecules()


if __name__ == "__main__":
    test_classifier()
