from catkit.surface import SlabGenerator
from ase.build import bulk


def test_surface_generator():
    atoms = bulk('Pd', 'fcc', a=4, cubic=True)

    # A variety of slab tests
    gen = SlabGenerator(
        atoms, miller_index=(1, 0, 0), layers=6, vacuum=4)

    gen = SlabGenerator(
        atoms, miller_index=(1, 1, 0), layers=6, vacuum=4)

    gen = SlabGenerator(
        atoms, miller_index=(0, 0, 1), layers=6, vacuum=4)

    # Test automatic termination finding
    atoms = bulk('Pd', 'fcc', a=4, cubic=True)
    atoms[3].symbol = 'Cu'

    gen = SlabGenerator(
        atoms, miller_index=(2, 1, 1), min_width=10, vacuum=4)

    gen.get_slab(iterm=1)


if __name__ == "__main__":
    test_surface_generator()
