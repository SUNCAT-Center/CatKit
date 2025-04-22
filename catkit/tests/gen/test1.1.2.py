from catkit.build import surface
from ase.io import write

# Size is 3 layers and 3 times the unit area of the primitve cell
atoms = surface('Pd', size=(3, 3), miller=(1, 1, 1), vacuum=4)

print(atoms.connectivity)

for i in atoms.get_surface_atoms():
    atoms[i].symbol = 'Au'

atoms.edit()
