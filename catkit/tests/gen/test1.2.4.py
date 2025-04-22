from catkit.gen.surface import SlabGenerator
from catkit.gen.adsorption import Builder
from ase.build import bulk
import numpy as np

atoms = bulk('Pd', 'fcc', a=4, cubic=True)
atoms[3].symbol = 'Cu'

gen = SlabGenerator(
    atoms,
    miller_index=(1, 1, 1),
    layers=4,
    fixed=2,
    vacuum=10)

slab = gen.get_slab()

builder = Builder(slab)
print(builder)
