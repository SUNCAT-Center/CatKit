import matplotlib
matplotlib.use('Agg')
from ase.io import read
from catgen.analysis import Classifier


def test_classifier():
    atoms = read('data/molecules.db')

    # This needs to be made rigorous, but also to be converted to
    # Gratoms automatically.
    cl = Classifier(atoms)
    ads_atoms = cl.id_adsorbate_atoms()
    slab_atoms = cl.id_slab_atoms()


if __name__ == "__main__":
    test_classifier()
