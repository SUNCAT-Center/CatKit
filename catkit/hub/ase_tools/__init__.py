from ase import Atoms
from ase.io import read
# from ase.io.trajectory import convert
import numpy as np
import ase
import copy
from catkit.hub.tools import get_atoms, get_state, clear_prefactor

# A lot of functions from os.path
# in python 2 moved to os. and changed
# their signature. Pathlib can be
# installed under python2.X with
# pip install pathlib2 and is in
# standard library in Python 3,
# hence we use it as a compatiblity
# library
try:
    from pathlib import Path
    Path().expanduser()
except (ImportError, AttributeError):
    from pathlib2 import Path

PUBLICATION_TEMPLATE = ''


def get_chemical_formula(atoms, mode='metal'):
    """
    Compatibility function, return mode=metal, when
    available, mode=hill, when not (ASE <= 3.13)
    """
    try:
        return atoms.get_chemical_formula(mode=mode)
    except ValueError:
        return atoms.get_chemical_formula(mode='hill')


def symbols(atoms):
    formula = get_chemical_formula(atoms)
    symbols = ase.atoms.string2symbols(formula)
    return ''.join(symbols)


def collect_structures(foldername, verbose=False, level='*'):
    structures = []
    if verbose:
        print(foldername)
    for i, filename in enumerate(Path(foldername).glob(level)):
        posix_filename = str(filename.as_posix())
        if verbose:
            print(i, posix_filename)
        if posix_filename.endswith('publication.txt'):
            with open(posix_filename) as infile:
                global PUBLICATION_TEMPLATE
                PUBLICATION_TEMPLATE = infile.read()
        elif Path(posix_filename).is_file():
            try:
                filetype = ase.io.formats.filetype(posix_filename)
            except Exception as e:
                continue
            if filetype:
                try:
                    structure = ase.io.read(posix_filename)
                    structure.info['filename'] = posix_filename
                    structure.info['filetype'] = ase.io.formats.filetype(
                        posix_filename)
                    try:
                        structure.get_potential_energy()
                        # ensure that the structure has an energy
                        structures.append(structure)
                    except RuntimeError:
                        print("Did not add {posix_filename} since it has no energy"
                              .format(
                                  posix_filename=posix_filename,
                              ))
                except StopIteration:
                    print("Warning: StopIteration {posix_filename} hit."
                          .format(
                              posix_filename=posix_filename,
                          ))
                except IndexError:
                    print("Warning: File {posix_filename} looks incomplete"
                          .format(
                              posix_filename=posix_filename,
                          ))
                except OSError as e:
                    print("Error with {posix_filename}: {e}".format(
                        posix_filename=posix_filename,
                        e=e,
                    ))
                except AssertionError as e:
                    print("Hit an assertion error with {posix_filename}: {e}".format(
                        posix_filename=posix_filename,
                        e=e,
                    ))
                except ValueError as e:
                    print("Trouble reading {posix_filename}: {e}".format(
                        posix_filename=posix_filename,
                        e=e,
                    ))
                except DeprecationWarning as e:
                    print("Trouble reading {posix_filename}: {e}".format(
                        posix_filename=posix_filename,
                        e=e,
                    ))
    return structures


def get_energies(atoms_list):
    """ Potential energy for a list of atoms objects"""
    if len(atoms_list) == 1:
        return atoms_list[0].get_potential_energy()
    elif len(atoms_list) > 1:
        energies = []
        for atoms in atoms_list:
            energies.append(atoms.get_potential_energy())
        return energies


def get_atomic_numbers(atoms):
    return list(atoms.get_atomic_numbers())


def get_formula_from_numbers(numbers):
    formula = Atoms(numbers).get_chemical_formula(mode='all')
    return formula


def get_numbers_from_formula(formula):
    atoms = Atoms(formula)
    return get_atomic_numbers(atoms)



def get_reaction_energy(structures, reaction, reaction_atoms, states,
                        prefactors, prefactors_TS, energy_corrections):
    energies = {}
    for key in structures.keys():
        energies.update({key: ['' for n in range(len(structures[key]))]})
    for key, atoms_list in structures.items():
        for i, atoms in enumerate(atoms_list):
            try:
                name = clear_prefactor(reaction[key][i])
            except BaseException:
                name = None
            if name in energy_corrections.keys():
                Ecor = energy_corrections[name]
            else:
                Ecor = 0
            energies[key][i] = prefactors[key][i] * \
                (atoms.get_potential_energy() + Ecor)

    # Reaction energy:
    energy_reactants = np.sum(energies['reactants'])
    energy_products = np.sum(energies['products'])

    reaction_energy = energy_products - energy_reactants

    # Activation energy
    if 'TS' in structures.keys():
        # Is a different empty surface used for the TS?
        if 'TSempty' in structure.keys():
            for key in reaction_atoms.keys():
                if '' in reaction_atoms[key]:
                    index = reaction_atoms[key].index('')
                    empty = structures[key][index]
            tsempty = structures['TSempty'][0]
            tsemptydiff = tsempty.get_potential_energy - \
                empty.get_potential_energy()

        for i, structure in enumerate(structures['reactants']):
            try:
                name = clear_prefactor(reaction['reactants'][i])
            except BaseException:
                name = None
            if name in energy_corrections.keys():
                Ecor = energy_corrections[name]
            else:
                Ecor = 0
            energies['reactants'][i] = prefactors_TS['reactants'][i]\
                * structure.get_potential_energy() + Ecor
            if 'TSempty' in structures.keys() and \
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
    tolerance = 0.01
    d = atoms.positions[:, 2]
    keys = np.argsort(d)
    ikeys = np.argsort(keys)
    mask = np.concatenate(([True], np.diff(d[keys]) > tolerance))
    layer_i = np.cumsum(mask)[ikeys]

    if layer_i.min() == 1:
        layer_i -= 1
    return layer_i


def get_surface_composition(atoms):
    if len(np.unique(atoms.get_atomic_numbers())) == 1:
        return atoms.get_chemical_symbols()[0]

    layer_i = get_layers(atoms)
    top_layer_i = np.max(layer_i)
    atom_i = np.where(layer_i >= top_layer_i - 1)[0]

    layer_atoms = atoms[atom_i]

    surface_composition = layer_atoms.get_chemical_formula(mode='metal')

    return surface_composition


def get_n_layers(atoms):
    layer_i = get_layers(atoms)
    n = np.max(layer_i)
    return n


def get_bulk_composition(atoms):
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


def check_in_ase(atoms, ase_db, energy=None):
    """Check if entry is allready in ASE db"""

    db_ase = ase.db.connect(ase_db)
    if energy is None:
        energy = atoms.get_potential_energy()
    formula = get_chemical_formula(atoms)
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


def write_ase(atoms, db_file, stdout, user=None, data=None, **key_value_pairs):
    """Connect to ASE db"""
    atoms = tag_atoms(atoms)
    db_ase = ase.db.connect(db_file)
    _normalize_key_value_pairs_inplace(key_value_pairs)
    id = db_ase.write(atoms, data=data, **key_value_pairs)
    stdout.write('  writing atoms to ASE db row id = {}\n'.format(id))
    unique_id = db_ase.get(id)['unique_id']
    return unique_id


def update_ase(db_file, identity, stdout, **key_value_pairs):
    """Connect to ASE db"""
    db_ase = ase.db.connect(db_file)

    _normalize_key_value_pairs_inplace(key_value_pairs)
    count = db_ase.update(identity, **key_value_pairs)
    stdout.write('  Updating {0} key value pairs in ASE db row id = {1}\n'
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
                sites.update({clear_prefactor(mol): site})
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
