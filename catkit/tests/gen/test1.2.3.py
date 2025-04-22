from catkit.gen.adsorption import AdsorptionSites
from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase import Atom
import numpy as np

bulk = bulk('Pd', 'fcc', a=5, cubic=True)
bulk[3].symbol = 'Cu'

gen = SlabGenerator(
    bulk,
    miller_index=(2, 1, 1),
    layers=10,
    vacuum=5)

atoms = gen.get_slab()
sites = AdsorptionSites(atoms)

coordinates = sites.get_coordinates()
vectors = sites.get_adsorption_vectors()

heights = np.arange(0, 2, 0.25)
for i, c in enumerate(coordinates):
    for h in heights:
        atoms += Atom('X', c + vectors[i] * h)

atoms.wrap()
atoms.edit()
