from catkit.gen.adsorption import AdsorptionSites
from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase.io import write
from ase import Atom

bulk = bulk('Pd', 'fcc', a=5, cubic=True)
bulk[3].symbol = 'Cu'

gen = SlabGenerator(
    bulk,
    miller_index=(3, 2, 1),
    layers=13,
    vacuum=5)

atoms = gen.get_slab(size=1)

sites = AdsorptionSites(atoms)

# Positon of each site
coordinates = sites.get_coordinates()

# Number of adjacent surface atoms
connectivity = sites.get_connectivity()

# The indices of adjacent surface atoms
topology = sites.get_topology()

# Only print every 5th entry.
print('Coordinates:\n', coordinates[::5], '\n')
print('Connectivity:\n', connectivity[::5], '\n')
print('Topology:\n', topology[::5], '\n')

periodic = sites.get_periodic_sites()
print('Sites by periodicity:\n', periodic[::5], '\n')

symmetric = sites.get_symmetric_sites()
print('Sites by symmetry:\n', symmetric[::5])

atm = {1: 'X', 2: 'He', 3: 'F', 4: 'N'}
for i, c in enumerate(coordinates):
    typ = connectivity[i]
    atoms += Atom(atm[typ], c + [0, 0, 2])

atoms.edit()
