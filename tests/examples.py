from catkit.surface import SlabGenerator
from ase.io import write
from ase.build import bulk
import unittest


class CatappBackendTestCase(unittest.TestCase):
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

        # Create side and top-down visuals
        img_name = './images/CuPd3-term-{}.pov'.format(i)
        write(
            img_name,
            slab,
            show_unit_cell=2,
            rotation='-90x',
        )

        write(
            img_name.replace('.pov', '-top.pov'),
            slab,
            show_unit_cell=2,
        )
