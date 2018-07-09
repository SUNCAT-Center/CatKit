from .fwio import array_to_list
import os
import contextlib
import json
import re
import hashlib
import msgpack
import numpy as np
from ase.io import read
from ase.db import connect


def geometry_hash(atoms):
    """A hash based strictly on the geometry features of
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
    """Does path management: if the path doesn't exists, create it
    otherwise, move into it until the indentation is broken.

    e.g.

    with cd('the/path/is/real'):
        'do things in the new path'

    Parameters
    ----------
    path : str
        Directory path to create and change into.
    """
    cwd = os.getcwd()
    try:
        if not os.path.exists(path):
            os.makedirs(path)
        os.chdir(path)
        yield
    finally:
        os.chdir(cwd)
