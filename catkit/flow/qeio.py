from .fwio import array_to_list
import os
import contextlib
import json
import re
import hashlib
import msgpack
import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator as SPC
from ase.db import connect
from ase.units import Ry, Bohr
from ase.io import read, write
# 1.889726 is the atomic unit of length per Angstrom (aul/A)
aul = 1.889726


def geometry_hash(atoms):
    """ A hash based strictly on the geometry features of
    an atoms object: positions, cell, and symbols.

    This is intended for planewave basis set calculations,
    so pbc is not considered.

    Each element is sorted in the algorithem to help prevent
    new hashs for identical geometries.
    """

    atoms = atoms.copy()
    atoms.wrap()

    pos = atoms.get_positions()

    # Sort the cell array by magnitude of z, y, x coordinates, in that order
    cell = np.array(sorted(atoms.get_cell(),
                           key=lambda x: (x[2], x[1], x[0])))

    # Flatten the array and return a string of numbers only
    # We only consider position changes up to 3 decimal places
    cell_hash = np.array_str(np.ndarray.flatten(cell.round(3)))
    cell_hash = ''.join(cell_hash.strip('[]').split()).replace('.', '')

    # Sort the atoms positions similarly, but store the sorting order
    pos = atoms.get_positions()
    srt = [i for i, _ in sorted(enumerate(pos),
                                key=lambda x: (x[1][2], x[1][1], x[1][0]))]
    pos_hash = np.array_str(np.ndarray.flatten(pos[srt].round(3)))
    pos_hash = ''.join(pos_hash.strip('[]').split()).replace('.', '')

    # Create a symbols hash in the same fashion conserving position sort order
    sym = np.array(atoms.get_atomic_numbers())[srt]
    sym_hash = np.array_str(np.ndarray.flatten(sym))
    sym_hash = ''.join(sym_hash.strip('[]').split())

    # Assemble a master hash and convert it through an md5
    master_hash = cell_hash + pos_hash + sym_hash
    md5 = hashlib.md5(master_hash.encode('utf-8'))
    _hash = md5.hexdigest()

    return _hash


def write_to_db(path, db_name='master.db', keys={}, traj=False, pdos=False):
    """ This function is used to write an ASE database entry for
    each atoms image inside a QE directory.
    """

    db = connect(os.path.join(os.getcwd(), db_name))

    # Change to the directory
    with cd(path):

        # Get the atoms object from the current directory
        if traj:
            images = read('output.traj', ':')
        else:
            images = [read('output.traj')]

        for i, atoms in enumerate(images):

            # Get keys-value-pairs from directory name.
            path = [x for x in os.getcwd().split('/') if '=' in x]

            keys.update({
                'traj': i,
                'rtraj': len(images) - i - 1})
            for key_value in path:
                key = key_value.split('=')[0]
                value = key_value.split('=')[1]

                # Try to recognize characters and convert to
                # specific data types for easy access later.
                if '.' in value:
                    value = float(value)
                elif value.isdigit():
                    value = int(value)
                elif value == 'False':
                    value = bool(False)
                elif value == 'True':
                    value = bool(True)
                else:
                    value = str(value)

                # Add directory keys
                keys[key] = value

            # Collect the psppath from the pw.inp file
            with open('pw.inp', 'r') as f:
                line = f.readline()
                while 'pseudo_dir' not in line:
                    line = f.readline()

                psppath = re.split("[']", line)[1]

            atoms.info['psppath'] = psppath

            data = {
                'path': os.getcwd(),
                'parameters': atoms.info,
                'geometry_hash': geometry_hash(atoms)}

            if pdos and i == 0:
                if os.path.exists('dos.msg'):
                    with open('dos.msg') as f:
                        dos_decrypt = msgpack.load(f)

                    dos_decrypt = list(dos_decrypt)
                    array_to_list(dos_decrypt)
                    data['dos'] = json.dumps(dos_decrypt, encoding='utf-8')

            # Generate/add-to the db file
            db.write(atoms, key_value_pairs=keys, data=data)


@contextlib.contextmanager
def cd(path):
    """ Does path management: if the path doesn't exists, create it
    otherwise, move into it.

    This is a context manager function, so it should be used as a with:
    e.x.

    with cd('the/path/is/real'):
        'do things in the new path'
    """

    cwd = os.getcwd()
    try:
        if not os.path.exists(path):
            os.makedirs(path)
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)


def attach_results(f, atoms, write_file=True):
    """ Return the TS corrected energy for a scf instance
    in a log file and attach them to the given atoms obejct.

    Will also attach the forces and stress if applicable.
    """
    energy, forces, stress = None, None, None

    line = f.readline()
    while '!    total energy' not in line:
        line = f.readline()

    energy = float(line.split()[-2]) * Ry

    # Correct for non-zero temperature smearing
    for i in range(20):

        line = f.readline()
        if '     smearing contrib.' in line:
            energy -= 0.5 * float(line.split()[-2]) * Ry

        # Collect the forces on the atoms
        if 'Forces acting on atoms (Ry/au):' in line:
            for _ in range(4):
                line = f.readline()
                if 'atom' in line:
                    break

            forces = []
            for _ in range(len(atoms)):
                forces += [line.split()[-3:]]
                line = f.readline()

            forces = np.array(forces, dtype=float) / Ry * aul

            # If forces were located, attempt to find stress
            for i in range(10):
                line = f.readline()

                if 'total   stress' in line:

                    stress = []
                    for _ in range(3):
                        line = f.readline()
                        stress += [line.split()[-3:]]

                    stress = np.array(stress, dtype=float) / Ry * Bohr ** 3
                    break

    # attach the calculator
    calc = SPC(
        atoms=atoms,
        energy=energy,
        forces=forces,
        stress=stress)
    atoms.set_calculator(calc)

    return atoms


def log_to_atoms(log_file='log', ent=-1, out_file=None):
    """ Parse a QE log file for atoms trajectory and return a list
    of atoms objects representative of the relaxation path.

    NOTE: trajectory information is only returned for calculations
    run with BFGS internal to QE.
    """

    images = []
    with open(log_file) as f:
        line = f.readline()

        # Flag to read trajectory 'ent' only
        with os.popen(
                'grep -n Giannozzi ' +
                log_file +
                ' 2>/dev/null', 'r') as p:
            n = int(p.readlines()[ent].split()[0].strip(':'))

        for i in range(n):
            line = f.readline()

        # Read lines one at a time
        while line:
            line = f.readline()

            # Signifies a new trajectory
            # Clear any existing values from previous runs
            if '(npk)' in line:

                # Look for an input trajectory in the same file and use it
                # Convenient for conserving constraints, tags, and atoms info
                in_file = os.path.join(
                    '/'.join(log_file.split('/')[:-1]),
                    'input.traj')

                if os.path.exists(in_file):

                    # If it does exist, read it in as the initial configuration
                    atoms = read(in_file)
                    atoms.wrap()
                    natoms = len(atoms)
                    pos = atoms.get_positions()

                    # Skip past the geometry information
                    while 'site n.' not in line:
                        line = f.readline()

                # Otherwise, collect from the data
                else:
                    atoms = None

            # Example properties
            ######################
            # bravais-lattice index     =            0
            # lattice parameter (alat)  =       1.8897  a.u.
            # unit-cell volume          =    3209.1777 (a.u.)^3
            # number of atoms/cell      =            5
            # number of atomic types    =            2
            # number of electrons       =        45.00
            # number of Kohn-Sham states=           55
            # kinetic-energy cutoff     =      36.7493  Ry
            # charge density cutoff     =     367.4932  Ry
            # convergence threshold     =      7.3E-08
            # mixing beta               =       0.1000
            # number of iterations used =            8  plain     mixing
            # Exchange-correlation      = BEEF ( 1  4 27 13 2)
            # nstep                     =           50

            # Collect potentially relevent properties
            # The elif can be omitted if order is assured
            elif 'number of atoms/cell      =' in line:
                natoms = int(line.split()[-1])

            # Collect cell dimensions
            elif 'celldm(1)' in line:
                alat = float(line.split()[1]) / aul

            elif 'crystal' in line:
                cell = []
                for _ in range(3):
                    line = f.readline()
                    cell += [[float(x) for x in line.split()[3:6]]]
                cell = np.array(cell) * alat

            # Collect positions, symbols, and number of atoms
            elif 'site n.' in line:
                pos, syms = [], []

                for _ in range(natoms):
                    line = f.readline()
                    pos += [line.split()[-4:-1]]
                    syms += [line.split()[1].strip('0123456789')]

                pos = np.array(pos, dtype=float) * alat

                # Setup the atoms object
                atoms = Atoms(syms, pos, cell=cell, pbc=(1, 1, 1))

            # This should be the last piece of information
            elif 'number of k points=' in line:

                atoms = attach_results(f, atoms)

                # Add atom to images
                images = [atoms]

                # Only atomic positions and energies need to be collected now
                # until the calculation ends
                while 'JOB DONE.' not in line and line:
                    line = f.readline()

                    # A duplicate of the coordinates printed previously
                    if 'Begin final coordinates' in line:
                        break

                    if 'ATOMIC_POSITIONS' in line:
                        atoms = atoms.copy()

                        coord = line.split('(')[-1]
                        for i in range(natoms):
                            line = f.readline()
                            pos[i][:] = line.split()[1:4]

                            # It's possible to recover constraints here,
                            # but not yet implemented If there are 7
                            # characters in the line, we have constraints
                            # if len(line.split()) == 7:
                            #     cons += [line.split()[-3:]]
                            # else:
                            #     cons += [[1] * 3]

                        # cons = np.array(cons, dtype=float)

                        if coord == 'alat)':
                            atoms.set_positions(pos * alat)
                        elif coord == 'bohr)':
                            atoms.set_positions(pos * Bohr)
                        elif coord == 'angstrom)':
                            atoms.set_positions(pos)
                        else:
                            atoms.set_scaled_positions(pos)

                        atoms = attach_results(f, atoms)
                        images += [atoms]

                if out_file:
                    write(out_file, images)

                return images
