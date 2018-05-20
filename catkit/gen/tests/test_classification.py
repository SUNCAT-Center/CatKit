from catkit.gen.analysis import Classifier
from ase.io import read
import pkg_resources


def test_classifier():
    data_path = pkg_resources.resource_filename(
        'catkit', 'data/molecules.db')
    atoms = read(data_path)

    # This needs to be made rigorous, but also to be converted to
    # Gratoms automatically.
    cl = Classifier(atoms)
    ads_atoms = cl.id_adsorbate_atoms()
    slab_atoms = cl.id_slab_atoms()


if __name__ == "__main__":
    test_classifier()
