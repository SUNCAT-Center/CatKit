from .fwio import array_to_list, atoms_to_encode
from espresso import Espresso
import ase
import msgpack
import json


def get_relaxed_calculation(in_file='output.traj'):
    """Attach a stored calculator in the current directory to the
    provided atoms object.

    Then return the atoms object with the calculator attached.

    Parameters
    ----------
    in_file : str
        Name of the relaxed trajectroy file to load.
    """
    atoms = ase.io.read(in_file)

    # Reinitialize the calculator from calc.tar.gz and attach it.
    calc = Espresso(**atoms.info)
    calc.load_flev_output()
    atoms.set_calculator(calc)

    return atoms


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


def get_pdos(out_file='dos.msg', resolution=2):
    """Calculate and save the projected DOS. Requires a previously relaxed
    calculation.

    Parameters
    ----------
    out_file : str
        Name of the output file to save the results to.
    resolution : int | array (3,)
        Multiple of the kpoints to increase their resolution by.
    """
    atoms = get_relaxed_calculation()
    calc = atoms.get_calculator()

    # Calculate the pdos and write to disk
    dos = calc.calc_pdos(
        nscf=True,
        kpts=atoms.info['kpts'] * resolution,
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
    images = ase.io.read(in_file, ':')

    for atoms in images[1:-1]:
        calc = Espresso(**atoms.info)
        atoms.set_calculator(calc)

    neb = ase.neb.NEB(images)
    opt = ase.optimize.BFGS(neb, trajectory='output.traj', logfile=None)
    opt.run(fmax=atoms.info.get('fmax'))
    out_images = ase.io.read('output.traj', ':')

    # Save the calculator to the local disk for later use.
    try:
        calc.save_flev_output()
    except(RuntimeError):
        calc.save_output()

    return atoms_to_encode(out_images)
