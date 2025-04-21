from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase import Atom

bulk = bulk('Pd', 'fcc', a=5, cubic=True)
bulk[3].symbol = 'Cu'

gen = SlabGenerator(
    bulk,
    miller_index=(1, 1, 1),
    layers=6,
    vacuum=10)

atoms = gen.get_slab()
coordinates, connectivity = gen.adsorption_sites(atoms)

atm = {1: 'X', 2: 'He', 3: 'F'}
for i, c in enumerate(coordinates):
    typ = connectivity[i]
    atoms += Atom(atm[typ], c + [0, 0, 2])

atoms.edit()
