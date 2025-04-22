from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase.visualize import view

# Make a test slab
atoms = bulk('Pd', 'fcc', cubic=True)
atoms[3].symbol = 'Cu'

gen = SlabGenerator(
    atoms,
    miller_index=(2, 1, 1),
    layers=9,
    fixed=5,
    vacuum=4)

terminations = gen.get_unique_terminations()

images = []
for i, t in enumerate(terminations):
    images += [gen.get_slab(iterm=i)]

view(images)
