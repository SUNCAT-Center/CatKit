from . import fwio
from . import utils
from . import db
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

    if calc.__name__ == 'Espresso':
        # Patch for reading magmom of trajectory
        images = decaf.io.read(out_file, ':')
    else:
        images = ase.io.read(out_file, ':')

    # Moneky patch for constraints and pbc conservation.
    for image in images:
        image.constraints = atoms.constraints
        image._pbc = atoms.pbc

    return fwio.atoms_to_encode(images)


def catflow_relaxation(calculator):
    """Performs a relaxation of an arbitrary structure in a manor
    consistent with catflow database conventions.

    The type of strucutre will be automatically identified with
    periodic boundry conditions even if unused by the calculator.

    Parameters
    ----------
    calculator : str
        String representation of a calculator import.
        This currently only supports an Espresso calculator (ASE or decaf).
    """
    atoms = ase.io.read('input.traj')
    parameters = atoms.info['calculator_parameters']
    data = atoms.info.copy()

    # Setting up the calculator
    calculator = utils.str_to_class(calculator)
    calc = calculator(atoms, **parameters)

    # Perform the calculation and write trajectory from log.
    atoms.get_potential_energy()

    if np.all(atoms.pbc):
        # If bulk, we need to resymmeterize and perform a final scf
        tag, atoms = ckutil.aflow.get_prototype_tag(atoms, tol=1e-2)
        parameters['calculation'] = 'scf'
        calc = calculator(atoms, **parameters)
        atoms.get_potential_energy()

    if calc.__name__ == 'Espresso':
        # Patch for reading magmom of trajectory
        images = decaf.io.read('pw.pwo', ':')
        if np.all(atoms.pbc):
            images += [atoms]
    else:
        images = ase.io.read('pw.pwo', ':')
        if np.all(atoms.pbc):
            images += [atoms]

    # Moneky patch for constraints and pbc conservation.
    for image in images:
        image.constraints = atoms.constraints
        image._pbc = atoms.pbc

    return fwio.atoms_to_encode(images)


def upload_relaxation():
    """Uploads a completed relaxation calculation to the Catflow database.
    This will require CATFLOW_DB to be specified on each calculation node
    running Fireworks.
    """
    with Connect() as db:
        db.
