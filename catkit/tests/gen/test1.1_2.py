from catkit.build import surface
from ase.build import bulk

# Make a test slab
atoms = bulk('Pd', 'fcc', cubic=True)
atoms[3].symbol = 'Cu'

slab = surface(atoms,
               size=(1, 1, 9), miller=(2, 1, 1),
               vacuum=9, fixed=5)

slab.edit()
