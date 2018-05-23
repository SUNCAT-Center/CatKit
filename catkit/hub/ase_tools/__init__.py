from ase import Atoms
from ase.io import read
import numpy as np
import ase
import copy


def read_ase(filename):
    import six
    if isinstance(filename, six.string_types):
        atoms = read(filename)
    else:
        atoms = filename
    return atoms


def check_traj(filename, strict=True, verbose=True):
    from ase.io.trajectory import convert
    import math
    try:
        atoms = read_ase(filename)
        if verbose:
            print('traj file read!')
    except BaseException:
        try:
            convert(filename)
            if verbose:
                print('Converting to new ase format!')
            atoms = read_ase(filename)
        except BaseException:
            print('Could not read traj file: {}'.format(filename))
            return False

    try:
        atoms.get_potential_energy()
        assert not math.isnan(atoms.get_potential_energy()), 'Energy is NaN!'
    except BaseException:
        if strict:
            raise RuntimeError('No energy for .traj file: {}'.format(filename))
        else:
            print('No energy for .traj file: {}'.format(filename))
            return False
    return True


def get_reference(filename):
    atoms = read_ase(filename)
    energy = atoms.get_potential_energy()
    name = atoms.get_chemical_formula()
    return {name: str(energy)}


def get_pbc(filename):
    atoms = read_ase(filename)
    return atoms.get_pbc()


def get_traj_str(filename):
    from ase.db.row import AtomsRow
    from ase.io.jsonio import encode
    atoms = read_ase(filename)
    row = AtomsRow(atoms)
    dct = {}
    for key in row.__dict__:
        if key[0] == '_' or key in row._keys or key == 'id':
            continue
        dct[key] = row[key]
    constraints = row.get('constraints')
    if constraints:
        dct['constraints'] = constraints

    txt = ','.join('"{0}": {1}'.format(key, encode(dct[key]))
                   for key in sorted(dct.keys()))

    atoms_txt = '{{{0}}}'.format(txt)
    return atoms_txt


def get_chemical_formula(filename, mode='metal'):
    atoms = read_ase(filename)
    return atoms.get_chemical_formula(mode=mode)


def get_number_of_atoms(filename):
    atoms = read_ase(filename)
    return atoms.get_number_of_atoms()


def get_energy_diff(filename, filename_ref):
    atoms = read_ase(filename)
    reference = read_ase(filename_ref)
    return atoms.get_potential_energy() - reference.get_potential_energy()


def get_energies(filenames):
    if len(filenames) == 1:
        atoms = read_ase(filenames[0])
        return atoms.get_potential_energy()
    elif len(filenames) > 1:
        energies = []
        for filename in filenames:
            atoms = read_ase(filename)
            energies.append(atoms.get_potential_energy())
        return energies


def get_energy(filename):
    atoms = read_ase(filename)
    return atoms.get_potential_energy()


def get_atomic_numbers(filename):
    atoms = read_ase(filename)
    return list(atoms.get_atomic_numbers())


def get_formula_from_numbers(numbers):
    formula = Atoms(numbers).get_chemical_formula(mode='all')
    return formula


def get_numbers_from_formula(formula):
    atoms = Atoms(formula)
    return get_atomic_numbers(atoms)


def clear_state(name):
    name = name.replace('*', '').replace('(g)', '')
    name = name.replace('star', '').replace('gas', '')
    return name


def clear_prefactor(molecule):
    if molecule == '':
        return molecule
    if not molecule[0].isalpha():
        i = 0
        while not molecule[i].isalpha():
            i += 1
        molecule = molecule[i:]
    return molecule


def get_atoms(molecule):
    molecule = clear_state(molecule)
    if molecule == '':
        prefactor = 1
        return molecule, prefactor
    try:
        return '', float(molecule)
    except BaseException:
        pass
    if not molecule[0].isalpha():
        i = 0
        while not molecule[i].isalpha():
            i += 1
        prefactor = molecule[:i]
        if prefactor == '-':
            prefactor = -1
        prefactor = float(prefactor)
        molecule = molecule[i:]
    else:
        prefactor = 1

    temp = ''
    for k in range(len(molecule)):
        if molecule[k].isdigit():
            for j in range(int(molecule[k]) - 1):
                temp += molecule[k - 1]
        else:
            temp += molecule[k]

    molecule = ''.join(sorted(temp))

    return molecule, prefactor


def get_state(name):
    if '*' in name or 'star' in name:
        state = 'star'
    elif 'gas' in name:
        state = 'gas'
    else:
        state = 'star'
    return state


def get_reaction_energy(traj_files, reaction, reaction_atoms, states,
                        prefactors, prefactors_TS, energy_corrections):
    energies = {}
    for key in traj_files.keys():
        energies.update({key: ['' for n in range(len(traj_files[key]))]})
    for key, trajlist in traj_files.items():
        for i, traj in enumerate(trajlist):
            try:
                trajname = clear_prefactor(reaction[key][i])
            except BaseException:
                trajname = None
            if trajname in energy_corrections.keys():
                Ecor = energy_corrections[trajname]
            else:
                Ecor = 0
            energies[key][i] = prefactors[key][i] * (get_energy(traj) + Ecor)

    # Reaction energy:
    energy_reactants = np.sum(energies['reactants'])
    energy_products = np.sum(energies['products'])

    reaction_energy = energy_products - energy_reactants

    # Activation energy
    if 'TS' in traj_files.keys():
        # Is a different empty surface used for the TS?
        if 'TSempty' in traj_files.keys():
            for key in reaction_atoms.keys():
                if '' in reaction_atoms[key]:
                    index = reaction_atoms[key].index('')
                    traj_empty = traj_files[key][index]
            traj_tsempty = traj_files['TSempty'][0]
            tsemptydiff = get_energy(traj_tsempty) - get_energy(traj_empty)

        for i, traj in enumerate(traj_files['reactants']):
            try:
                trajname = clear_prefactor(reaction['reactants'][i])
            except BaseException:
                trajname = None
            if trajname in energy_corrections.keys():
                Ecor = energy_corrections[trajname]
            else:
                Ecor = 0
            energies['reactants'][i] = prefactors_TS['reactants'][i]\
                * (get_energy(traj) + Ecor)
            if 'TSempty' in traj_files.keys() and \
                    states['reactants'][i] == 'star':
                energies['reactants'][i] += prefactors_TS['reactants'][i]\
                    * tsemptydiff
        energy_reactants = np.sum(energies['reactants'])
        energy_TS = energies['TS'][0]
        activation_energy = energy_TS - energy_reactants
    else:
        activation_energy = None

    return reaction_energy, activation_energy


def tag_atoms(atoms, types=None):
    non_metals = ['H', 'He', 'B', 'C', 'N', 'O', 'F', 'Ne',
                  'Si', 'P', 'S', 'Cl', 'Ar',
                  'Ge', 'As', 'Se', 'Br', 'Kr',
                  'Sb', 'Te', 'I', 'Xe',
                  'Po', 'At', 'Rn']

    layer_i = get_layers(atoms)
    top_layer_i = np.max(layer_i)
    i = 0

    for i in range(0, top_layer_i + 1):
        atoms_i = np.where(layer_i == top_layer_i - i)[0]
        if len(np.where(layer_i == top_layer_i - i)[0]) == 1 and i < 4:
            atom = atoms[atoms_i[0]]
            if types is not None:
                if atom.symbol in types:
                    atom.tag = 0
            elif types is None:
                if atom.symbol in non_metals:
                    atom.tag = 0
        else:
            for l in atoms_i:
                atoms[l].tag = i + 1

    return atoms


def get_layers(atoms):
    # WARNING: this function is defined twice
    # with different parameter choices
    tolerance = 0.2
    d = atoms.positions[:, 2]
    keys = np.argsort(d)
    ikeys = np.argsort(keys)
    mask = np.concatenate(([True], np.diff(d[keys]) > tolerance))
    layer_i = np.cumsum(mask)[ikeys]

    if layer_i.min() == 1:
        layer_i -= 1
    return layer_i


def get_layers(atoms):
    # WARNING: this function is defined twice
    # with different parameter choices
    tolerance = 0.01
    d = atoms.positions[:, 2]
    keys = np.argsort(d)
    ikeys = np.argsort(keys)
    mask = np.concatenate(([True], np.diff(d[keys]) > tolerance))
    layer_i = np.cumsum(mask)[ikeys]

    if layer_i.min() == 1:
        layer_i -= 1
    return layer_i


def get_surface_composition(filename):
    filename = "{]".format(filename)
    atoms = read_ase(filename)

    if len(np.unique(atoms.get_atomic_numbers())) == 1:
        return atoms.get_chemical_symbols()[0]

    layer_i = get_layers(atoms)
    top_layer_i = np.max(layer_i)
    atom_i = np.where(layer_i >= top_layer_i - 1)[0]

    layer_atoms = atoms[atom_i]

    surface_composition = layer_atoms.get_chemical_formula(mode='metal')

    return surface_composition


def get_n_layers(filename):
    atoms = read_ase(filename)
    layer_i = get_layers(atoms)
    n = np.max(layer_i)
    return n


def get_bulk_composition(filename):
    atoms = read_ase(filename)

    if len(np.unique(atoms.get_atomic_numbers())) == 1:
        return atoms.get_chemical_symbols()[0]

    layer_i = get_layers(atoms)
    top_layer_i = np.max(layer_i)
    compositions = []
    for i in range(0, top_layer_i + 1):
        atom_i = np.where(layer_i == top_layer_i - i)[0]
        atoms_layer = atoms[atom_i]
        if len(np.unique(atoms_layer.get_atomic_numbers())) == 1:
            c = atoms_layer.get_chemical_symbols()[0]
            compositions.append(c)
        else:
            c = atoms[atom_i].get_chemical_formula(mode='metal')
            compositions.append(c)

    compositions = np.array(compositions)
    same_next_layer = compositions[1:] == compositions[:-1]
    bulk_compositions = compositions[:-1][same_next_layer]

    if len(bulk_compositions) > 0 and \
       all(c == bulk_compositions[0] for c in bulk_compositions):
        bulk_composition = bulk_compositions[0]
    else:
        bulk_composition = None
    return bulk_composition


def check_in_ase(filename, ase_db, energy=None):
    """Check if entry is allready in ASE db"""

    db_ase = ase.db.connect(ase_db)
    atoms = read_ase(filename)
    if energy is None:
        energy = atoms.get_potential_energy()
    formula = atoms.get_chemical_formula(mode='metal')
    rows = db_ase.select(energy=energy)
    n = 0
    ids = []
    for row in rows:
        if formula == row.formula:
            n += 1
            ids.append(row.id)
    if n > 0:
        id = ids[0]
        unique_id = db_ase.get(id)['unique_id']
        return id, unique_id
    else:
        return None, None


def _normalize_key_value_pairs_inplace(data):
    for key in data:
        if isinstance(data[key], np.int64):
            data[key] = int(data[key])


def write_ase(filename, db_file, user=None, data=None, **key_value_pairs):
    """Connect to ASE db"""
    atoms = read_ase(filename)
    atoms = tag_atoms(atoms)
    db_ase = ase.db.connect(db_file)
    _normalize_key_value_pairs_inplace(key_value_pairs)
    id = db_ase.write(atoms, data=data, **key_value_pairs)
    print('writing atoms to ASE db row id = {}'.format(id))
    unique_id = db_ase.get(id)['unique_id']
    return unique_id


def update_ase(db_file, identity, **key_value_pairs):
    """Connect to ASE db"""
    db_ase = ase.db.connect(db_file)

    _normalize_key_value_pairs_inplace(key_value_pairs)
    count = db_ase.update(identity, **key_value_pairs)
    print('Updating {0} key value pairs in ASE db row id = {1}'
          .format(count, identity))
    return


def get_reaction_from_folder(folder_name):
    reaction = {}
    if '__' in folder_name:  # Complicated reaction
        if '-' in folder_name and '_-' not in folder_name:
            # intermediate syntax
            a, b = folder_name.split('-')
            folder_name = a + '_-' + b

        reaction.update({'reactants': folder_name.split('__')[0].split('_'),
                         'products': folder_name.split('__')[1].split('_')})

    elif '_' in folder_name:  # Standard format
        AB, A, B = folder_name.split('_')
        if '-' in A:
            A = A.split('-')
            A[1] = '-' + A[1]
            products = [A[0], A[1], B]
        else:
            products = [A, B]
        reaction.update({'reactants': [AB],
                         'products': products})
    else:
        raise AssertionError('problem with folder {}'.format(folder_name))

    sites = {}
    for key, mollist in reaction.items():
        for n, mol in enumerate(mollist):
            if '@' in mol:
                mol, site = mol.split('@')
                sites.update({mol: site})
                reaction[key][n] = mol
            if 'gas' not in mol and 'star' not in mol:
                reaction[key][n] = mol + 'star'

    for key, mollist in reaction.items():
        n_star = mollist.count('star')
        if n_star > 1:
            for n in range(n_star):
                mollist.remove('star')
            mollist.append(str(n_star) + 'star')

    return reaction, sites


def get_reaction_atoms(reaction):
    reaction_atoms = {'reactants': [],
                      'products': []}

    prefactors = {'reactants': [],
                  'products': []}

    states = {'reactants': [],
              'products': []}

    for key, mollist in reaction.items():
        for molecule in mollist:
            atoms, prefactor = get_atoms(molecule)
            reaction_atoms[key].append(atoms)
            prefactors[key].append(prefactor)
            state = get_state(molecule)
            states[key].append(state)

    prefactors_TS = copy.deepcopy(prefactors)

    # Balance the number of slabs on each side of reaction
    n_star = {'reactants': 0,
              'products': 0}

    for key, statelist in states.items():
        for j, s in enumerate(statelist):
            if s == 'star':
                n_star[key] += prefactors[key][j]

    n_r = n_star['reactants']
    n_p = n_star['products']

    diff = n_p - n_r
    if abs(diff) > 0:
        if diff > 0:  # add empty slabs to left-hand side
            n_r += diff
            key = 'reactants'
        else:  # add to right-hand side
            diff *= -1  # diff should be positive
            n_p += diff
            key = 'products'

        if '' not in reaction_atoms[key]:
            reaction[key].append('star')
            prefactors[key].append(diff)
            if key == 'reactants':
                prefactors_TS[key].append(1)
            states[key].append('star')
            reaction_atoms[key].append('')
        else:
            index = states[key].index('star')
            prefactors[key][index] += diff
            if key == 'reactants':
                prefactors_TS[key][index] += diff

    if n_r > 1:  # Balance slabs for transition state
        count_empty = 0
        if '' in reaction_atoms['reactants']:
            index = reaction_atoms['reactants'].index('')
            count_empty = prefactors_TS['reactants'][index]
            prefactors_TS['reactants'][index] = -(n_r - count_empty - 1)
        else:
            reaction_atoms['reactants'].append('')
            prefactors['reactants'].append(0)
            states['reactants'].append('star')
            prefactors_TS['reactants'].append(-n_r + 1)
    else:
        if '' in reaction_atoms['reactants']:
            index = reaction_atoms['reactants'].index('')
            prefactors_TS['reactants'][index] = 1

    return reaction_atoms, prefactors, prefactors_TS, states


def debug_assert(expression, message, debug=False):
    if debug:
        try:
            assert expression, message
        except AssertionError as e:
            print(e)
            return False
    else:
        assert expression, message

    return True
