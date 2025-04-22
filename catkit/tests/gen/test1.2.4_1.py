from catkit.build import molecule
from catkit.gen.surface import SlabGenerator
from catkit.gen.adsorption import Builder
from ase.build import bulk

atoms = bulk('Pd', 'fcc', a=4, cubic=True)
atoms[3].symbol = 'Cu'

gen = SlabGenerator(
    atoms,
    miller_index=(1, 1, 1),
    layers=4,
    vacuum=4)

slab = gen.get_slab()

adsorbate = molecule('CH3')[0]
adsorbate.set_tags([-1, 0, 0, 0])

builder = Builder(slab)
ads_slab = builder.add_adsorbate(adsorbate, index=0)

ads_slab.edit()
