from catkit.gen.surface import SlabGenerator
from ase.build import bulk
from ase.visualize import view
from ase.io import read
from sys import argv


# First run get_restart.py in bulk dir to make sure restart.json is optimized structure

atoms = read('restart.json')

N = 8

gen = SlabGenerator(
    atoms,
    miller_index=(1, 1, 0),
    layers= N,
    layer_type='trim',
    standardize_bulk=True,
    symmetric=True,
    stoich=True,
    fixed=0,
    vacuum=7.5)


slabs = gen.get_slabs()

view(slabs)

