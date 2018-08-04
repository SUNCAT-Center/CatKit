from .fwio import atoms_to_encode
import functools
import ase
import sys
try:
    import espresso
except(ImportError):
    pass


def str_to_class(classname):
    return functools.reduce(
        getattr, classname.split('.'), sys.modules[__name__])


def get_potential_energy(
        calculator,
        in_file='input.traj',
        out_file='pw.pwo'):
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
    out_file : str
        Name of the output file to read the completed trajectory form.
    """
    atoms = ase.io.read(in_file)

    # Planewave basis set requires periodic boundary conditions
    atoms.set_pbc([1, 1, 1])

    # Setting up the calculator
    calculator = str_to_class(calculator)
    calculator(atoms, **atoms.info)

    # Perform the calculation and write trajectory from log.
    atoms.get_potential_energy()

    images = ase.io.read(out_file, ':')
    for image in images:
        image.constraints = atoms.constraints

    return atoms_to_encode(images)
