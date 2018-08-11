from . import fwio
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

    # Setting up the calculator
    calculator = str_to_class(calculator)
    calc = calculator(atoms, **atoms.info)

    # Perform the calculation and write trajectory from log.
    atoms.get_potential_energy()

    if isinstance(calc, espresso.Espresso):
        # Patch for reading magmom of trajectory
        images = espresso.io.read(out_file, ':')
        if 'magmoms' in calc.results:
            calc.get_pdos(update_projections=True)
            images[-1]._calc.results = calc.results.copy()
    else:
        images = ase.io.read(out_file, ':')

    # Moneky patch for constraints
    for image in images:
        image.constraints = atoms.constraints
        image._pbc = atoms.pbc

    return fwio.atoms_to_encode(images)
