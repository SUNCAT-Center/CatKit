from ase.calculators.singlepoint import SinglePointCalculator as SPC
from ase.constraints import dict2constraint
from ase.io import write
from ase import Atoms
import numpy as np
import json
supported_properties = ['energy', 'forces', 'stress', 'magmoms', 'magmom']


def array_to_list(data):
    """A function to covert all arrays in a structure of
    embeded dictionaries and lists into lists themselves.
    """
    if isinstance(data, list):
        for i, v in enumerate(data):
            print(i, v)
            if isinstance(v, np.ndarray):
                data[i] = v.tolist()

            elif isinstance(v, dict):
                array_to_list(v)
            elif isinstance(v, list):
                array_to_list(v)

    elif isinstance(data, dict):
        for k, v in list(data.items()):
            if isinstance(v, np.ndarray):
                data[k] = v.tolist()

            elif isinstance(v, dict):
                array_to_list(v)
            elif isinstance(v, list):
                array_to_list(v)


def encode_to_atoms(encode, out_file='input.traj'):
    """Dump the encoding to a local traj file."""
    # First, decode the trajectory
    data = json.loads(encode, encoding='utf-8')

    # Construct the initial atoms object
    atoms = Atoms(
        data['numbers'],
        data['trajectory']['0']['positions'],
        cell=data['trajectory']['0']['cell'],
        pbc=data['pbc'])
    atoms.info['calculator_parameters'] = data['calculator_parameters']
    atoms.set_constraint([dict2constraint(_) for _ in data['constraints']])
    initial_magmoms = data.get('initial_magmoms')
    if initial_magmoms:
        atoms.set_initial_magnetic_moments(initial_magmoms)

    # Attach the calculator
    results = {'atoms': atoms}
    for prop in supported_properties:
        results.update({prop: data['trajectory']['0'].get(prop)})
    calc = SPC(**results)
    atoms.set_calculator(calc)

    # Collect the rest of the trajectory information
    images = [atoms]
    for i in range(len(data['trajectory']))[1:]:
        atoms = atoms.copy()

        if data['trajectory'][str(i)]['cell']:
            atoms.set_cell(data['trajectory'][str(i)]['cell'])

        if data['trajectory'][str(i)]['positions']:
            atoms.set_positions(data['trajectory'][str(i)]['positions'])

        results = {'atoms': atoms}
        for prop in supported_properties:
            results.update({prop: data['trajectory'][str(i)].get(prop)})
        calc = SPC(**results)
        atoms.set_calculator(calc)

        images += [atoms]

    # Write the traj file
    if out_file:
        write(out_file, images)

    return images


def atoms_to_encode(images):
    """Converts an list of atoms objects to an encoding
    from a .traj file.
    """
    if not isinstance(images, list):
        images = [images]

    # Convert all constraints into dictionary format
    constraints = [_.todict() for _ in images[0].constraints]
    for i, C in enumerate(constraints):

        # Turn any arrays in the kwargs into lists
        for k, v in list(C['kwargs'].items()):
            if isinstance(v, np.ndarray):
                constraints[i]['kwargs'][k] = v.tolist()

    # Convert any arrays from the parameter settings into lists
    keys = images[0].info['calculator_parameters']
    array_to_list(keys)

    data = {'trajectory': {}}
    # Assemble the compressed dictionary of results
    for i, atoms in enumerate(images):

        if i == 0:
            # For first images, collect cell and positions normally
            pos = atoms.get_positions()
            update_pos = pos

            cell = atoms.get_cell()
            update_cell = cell

            # Add the parameters which do not change
            data['numbers'] = images[0].get_atomic_numbers().tolist()
            data['pbc'] = images[0].get_pbc().tolist()
            data['constraints'] = constraints
            data['calculator_parameters'] = keys
            initial_magmoms = atoms.arrays.get('initial_magmoms')
            if initial_magmoms is not None:
                data['initial_magmoms'] = list(initial_magmoms)

        else:
            # For consecutive images, check for duplication
            # If duplicates are found, do not store it
            if np.array_equal(atoms.get_positions(), pos):
                update_pos = np.array([])
            else:
                pos = atoms.get_positions()
                update_pos = pos

            if np.array_equal(atoms.get_cell(), cell):
                update_cell = np.array([])
            else:
                cell = atoms.get_cell()
                update_cell = cell

        results = {'positions': update_pos, 'cell': update_cell}
        if atoms._calc:
            for prop in supported_properties:
                results.update({prop: atoms._calc.results.get(prop)})

        for k, v in results.items():
            if isinstance(v, np.ndarray):
                results[k] = v.tolist()

        # Store trajectory, throwing out None values
        data['trajectory'][i] = {
            k: v for k, v in list(
                results.items()) if v is not None}

    # Return the reduced results in JSON compression
    encoding = json.dumps(data)

    return encoding
