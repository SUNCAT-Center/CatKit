from . import fwio
from . import utils
from . import db
import numpy as np
import ase
try:
    import decaf
    import ckutil
except(ImportError):
    pass


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
    calculator = utils.str_to_class(calculator)
    calc = calculator(atoms, **atoms.info['calculator_parameters'])

    # Perform the calculation and write trajectory from log.
    atoms.get_potential_energy()

    if isinstance(calc, decaf.Espresso):
        # Patch for reading magmom of trajectory
        images = decaf.io.read(out_file, ':')
    else:
        images = ase.io.read(out_file, ':')

    # Moneky patch for constraints and pbc conservation.
    images[0].info = atoms.info
    for image in images:
        image.constraints = atoms.constraints
        image._pbc = atoms.pbc

    return fwio.atoms_to_encode(images)


def catflow_relaxation(atoms=None, calculator_name=None, parameters=None):
    """Performs a relaxation of an arbitrary structure in a manor
    consistent with catflow database conventions.

    The type of strucutre will be automatically identified with
    periodic boundry conditions even if unused by the calculator.

    Parameters
    ----------
    calculator : str
        String representation of a calculator import.
        This currently only supports an Espresso calculator (ASE or decaf).
    parameters : dict
        Calculation parameters to use.
    """
    if atoms is None:
        atoms = ase.io.read('input.traj')
    if parameters is None:
        parameters = atoms.info['calculator_parameters']
    if calculator_name is None:
        calculator_name = parameters.pop('calculator_name')
    data = atoms.info.copy()

    # Setting up the calculator
    calculator = utils.str_to_class(calculator_name)
    calc = calculator(atoms, calc_file='pw0.pwi', **parameters)

    # Perform the calculation and write trajectory from log.
    atoms.get_potential_energy()

    if np.all(atoms.pbc):
        # If bulk, we need to resymmeterize and perform a final scf
        tag, atoms = ckutil.aflow.get_prototype_tag(calc.atoms, tol=1e-2)
        parameters['calculation'] = 'scf'
        calc = calculator(atoms, calc_file='pw1.pwi', **parameters)
        atoms.get_potential_energy()

    if isinstance(calc, decaf.Espresso):
        # Patch for reading magmom of trajectory
        images = decaf.io.read('pw0.pwo', ':')
        images[0].info = data
        if np.all(atoms.pbc):
            images += [atoms]
    else:
        images = ase.io.read('pw0.pwo', ':')
        images[0].info = data
        if np.all(atoms.pbc):
            images += [atoms]

    # Moneky patch for constraints and pbc conservation.
    for image in images:
        image.constraints = atoms.constraints
        image._pbc = atoms.pbc

    with db.Connect() as dbflow:
        dbflow.update_bulk_entry(images)

    return fwio.atoms_to_encode(images)


def get_bayesian_neb(
        calculator,
        neb_class,
        images='input.traj',
        out_file='neb.json'):
    """Performs a Bayesian NEB with a compatible ASE calculator.
    Keywords are defined inside the atoms object information.

    Parameters
    ----------
    calculator : str
        Import name of the calculator to use.
    neb_class : str
        Import name of the Bayesian optimizer to use.
    images: list of Atoms objects | str
        User provided iniital guess for the NEB pathway. Parameters for
        the calculator should be placed into the first image.
    out_file: str
        Name of the file to store output.
    """
    if isinstance(images, str):
        images = ase.io.read(images, ':')

    parameters = images[0].info['calculator_parameters']
    parameters['calculation'] = 'scf'
    parameters['tprnfor'] = True

    calculator = utils.str_to_class(calculator)
    neb_class = utils.str_to_class(neb_class)

    neb = neb_class(
        trajectory=images,
        calculator=calculator(**parameters),
    )

    fmax = parameters.get('fmax', 0.05)
    neb.run(fmax=fmax)

    images = ase.io.read(out_file, ':')
    images[0].info = parameters

    return fwio.atoms_to_encode(images)
