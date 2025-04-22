from catkit.build import molecule
from catkit.gen.surface import SlabGenerator
from catkit.gen.adsorption import Builder
from ase.build import bulk
from ase.visualize import view

atoms = bulk('Pd', 'fcc', a=4, cubic=True)
atoms[3].symbol = 'Cu'

gen = SlabGenerator(
    atoms,
    miller_index=(1, 1, 1),
    layers=4,
    vacuum=4)

slab = gen.get_slab()
adsorbate = molecule('C2H3')[1]

builder = Builder(slab)
ads_slab = builder.add_adsorbate(adsorbate, bonds=[0, 1], index=-1)

print('{} adsorption structures generated'.format(len(ads_slab)))

view(ads_slab)
