from .fwio import array_to_list, atoms_to_encode
from .qeio import log_to_atoms
from .hpcio import get_nnodes
from ase.dft.bee import BEEFEnsemble
from ase.io import read, write
from ase.neb import NEB
from ase.optimize import BFGS
import msgpack
import numpy as np
import json
try:
    from espresso import espresso
except(ImportError):
    pass


def get_relaxed_calculation(in_file='output.traj'):
    """Attach a stored calculator in the current directory to the
    provided atoms object.

    Then return the atoms object with the calculator attached.

    Parameters
    ----------
    in_file : str
        Name of the relaxed trajectroy file to load.
    """
    atoms = read(in_file)

    # Reinitialize the calculator from calc.tgz and attach it.
    calc = espresso(**atoms.info)
    calc.load_flev_output()
    atoms.set_calculator(calc)

    return atoms


def get_potential_energy(in_file='input.traj'):
    """Performs a ASE get_potential_energy() call with the ase-espresso calculator
    with the keywords defined inside the atoms object information.

    This can be a singlepoint calculation or a full relaxation depending
    on the keywords.

    Parameters
    ----------
    in_file : str
        Name of the input file to load from the local directory.
    """
    atoms = read(in_file)

    # Planewave basis set requires periodic boundary conditions
    atoms.set_pbc([1, 1, 1])

    # Assign kpoints to be split across nodes
    if get_nnodes() > 1:
        if not np.prod(atoms.info['kpts']) == 1:
            atoms.info['parflags'] = '-npool {}'.format(get_nnodes())

    # Setting up the calculator
    calc = espresso(**atoms.info)
    atoms.set_calculator(calc)

    # Perform the calculation and write trajectory from log.
    atoms.get_potential_energy()
    images = log_to_atoms(out_file='output.traj')

    # Save the calculator to the local disk for later use.
    try:
        calc.save_flev_output()
    except(RuntimeError):
        calc.save_output()

    if images[-1].info.get('beefensemble'):
        beef = BEEFEnsemble(calc).get_ensemble_energies()
        images[-1].info['beef_std'] = beef.std()
        write('output.traj', images)

    return atoms_to_encode(images)


def get_total_potential(out_file='potential.msg'):
    """Calculate and save the total potential. Requires a previously relaxed
    calculation.

    Parameters
    ----------
    out_file : str
        Name of the output file to save the results to.
    """
    atoms = get_relaxed_calculation()
    calc = atoms.get_calculator()

    # Collect the total potential and write to disk
    potential = calc.extract_total_potential()

    potential = list(potential)
    array_to_list(potential)

    # If outfile, write a MessagePack encoded version to disk
    if out_file:
        with open(out_file, 'w') as f:
            msgpack.dump(potential, f)

    # Return a BSON friendly version
    return json.dumps(potential, encoding='utf-8')


def get_pdos(out_file='dos.msg'):
    """Calculate and save the projected DOS. Requires a previously relaxed
    calculation.

    Parameters
    ----------
    out_file : str
        Name of the output file to save the results to.
    """
    atoms = get_relaxed_calculation()
    calc = atoms.get_calculator()

    # Calculate the pdos and write to disk
    dos = calc.calc_pdos(
        nscf=True,
        kpts=atoms.info['kpts'] * [2, 2, 1],
        DeltaE=0.01,
        slab=True,
        Emin=-40,
        Emax=40,
        tetrahedra=False,
        sigma=0.2)

    dos = list(dos)
    array_to_list(dos)

    # If outfile, write a MessagePack encoded version to disk
    if out_file:
        with open(out_file, 'w') as f:
            msgpack.dump(dos, f)

    # Return a BSON friendly version
    return json.dumps(dos, encoding='utf-8')


def get_neb(in_file='input.traj'):
    """Performs a ASE NEB optimization with the ase-espresso
    calculator with the keywords defined inside the atoms object information.

    Parameters
    ----------
    in_file : str
        Name of the input file to load from the local directory.
    """
    images = read(in_file, ':')

    for atoms in images:
        atoms.set_pbc([1, 1, 1])
        calc = espresso(**atoms.info)
        atoms.set_calculator(calc)

    neb = NEB(images)
    opt = BFGS(neb, trajectory='output.traj', logfile=None)
    opt.run(fmax=0.05)
    out_images = read('output.traj')

    # Save the calculator to the local disk for later use.
    try:
        calc.save_flev_output()
    except(RuntimeError):
        calc.save_output()

    return atoms_to_encode(out_images)
