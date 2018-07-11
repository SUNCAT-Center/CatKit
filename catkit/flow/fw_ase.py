from .hpcio import get_nnodes
from .fwio import atoms_to_encode
import numpy as np
import functools
import espresso
import ase
import sys


def str_to_class(classname):
    return functools.reduce(
        getattr, classname.split('.'), sys.modules[__name__])


def get_potential_energy(calculator, in_file='input.traj'):
    """Performs a ASE get_potential_energy() call with a compatible ase calculator.
    Keywords are defined inside the atoms object information.

    This can be a singlepoint calculation or a full relaxation depending
    on the keywords.

    Parameters
    ----------
    calculator : str
        String representation of a calculator import
    in_file : str
        Name of the input file to load from the local directory.
    """
    atoms = ase.io.read(in_file)

    # Planewave basis set requires periodic boundary conditions
    atoms.set_pbc([1, 1, 1])

    # Assign kpoints to be split across nodes
    if get_nnodes() > 1:
        if not np.prod(atoms.info['kpts']) == 1:
            atoms.info['parflags'] = '-npool {}'.format(get_nnodes())

    # Setting up the calculator
    calculator = str_to_class(calculator)
    calculator(atoms, **atoms.info)

    # Perform the calculation and write trajectory from log.
    atoms.get_potential_energy()
    images = ase.io.read('pwcsf.pwo', ':')

    return atoms_to_encode(images)
