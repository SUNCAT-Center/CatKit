from .cathubsqlite import CathubSQLite
from .tools import get_bases
from .import ase_tools

from datetime import date
import numpy as np
import os
import copy
import json


class FolderReader:
    """
    Class for reading data from organized folders and writing to local
    CathubSQLite database. Folders should be arranged with
    make_folders_template and are read in the order:

    level:

    0    folder_name
    1    |-- publication
    2        |-- dft_code
    3            |-- dft_functional
    4                |-- gas
    4                |-- metal1
    5                    |-- facet
    6                        |-- reaction

    Parameters
    ----------
    foldername: str
    debug: bool
        default is False. Choose True if the folderreader should  continue
        in spite of errors.
    update: bool
        Update data if allready present in database file. defalt is True
    """

    def __init__(self, folder_name, debug=False, strict=True, verbose=False,
                 update=True):
        self.debug = debug
        self.strict = strict
        self.verbose = verbose
        self.update = update

        self.catbase, self.data_base, self.user, self.user_base \
            = get_bases(folder_name=folder_name)
        self.user_base_level = len(self.user_base.split("/"))

        self.pub_level = 1
        self.DFT_level = 2
        self.XC_level = 3
        self.reference_level = 4
        self.slab_level = 5
        self.reaction_level = 6
        self.final_level = 6

    def read(self, skip=[], goto_metal=None, goto_reaction=None):
        """
        Get reactions from folders.

        Parameters
        ----------
        skip: list of str
            list of folders not to read
        goto_reaction: str
            Skip ahead to this metal
        goto_reaction:
            Skip ahead to this reacion
        """
        self.omit_folders = []
        self.coverages = None
        if len(skip) > 0:
            for skip_f in skip:
                self.omit_folders.append(skip_f)

        found_reaction = False
        for root, dirs, files in os.walk(self.user_base):
            for omit_folder in self.omit_folders:  # user specified omit_folder
                if omit_folder in dirs:
                    dirs.remove(omit_folder)
            level = len(root.split("/")) - self.user_base_level

            if level == self.pub_level:
                self.read_pub(root)

            if level == self.DFT_level:
                self.DFT_code = read_name_from_folder(root)

            if level == self.XC_level:
                self.DFT_functional = read_name_from_folder(root)
                self.read_gas(root + '/gas/')

            if level == self.reference_level:
                if 'gas' in root.split("/")[-1]:
                    continue

                if goto_metal is not None:
                    if root.split("/")[-1] == goto_metal:
                        goto_metal = None
                    else:
                        dirs[:] = []  # don't read any sub_dirs
                        continue
                self.read_bulk(root, files)

            if level == self.slab_level:
                self.read_slab(root, files)

            if level == self.reaction_level:
                if goto_reaction is not None:
                    if root.split("/")[-1] == goto_reaction:
                        goto_reaction = None
                    else:
                        dirs[:] = []  # don't read any sub_dirs
                        continue

                self.read_reaction(root, files)

            if level == self.final_level:
                self.read_energies(root, files)
                if self.key_value_pairs_reaction is not None:
                    yield self.key_value_pairs_reaction

    def write(self, skip=[], goto_reaction=None):
        for key_values in self.read(skip=skip, goto_reaction=goto_reaction):
            with CathubSQLite(self.cathub_db) as db:
                id = db.check(
                    key_values['chemical_composition'],
                    key_values['reaction_energy'])
                if id is None:
                    id = db.write(key_values)
                    print('Written to reaction db row id = {}'.format(id))
                elif self.update:
                    db.update(id, key_values)
                    print('Updated reaction db row id = {}'.format(id))
                else:
                    print('Already in reaction db with row id = {}'.format(id))

    def write_publication(self, pub_data):
        with CathubSQLite(self.cathub_db) as db:
            pid = db.check_publication(self.pub_id)
            if pid is None:
                pid = db.write_publication(pub_data)
                print('Written to publications db row id = {}'.format(pid))
        return pid

    def read_pub(self, root):
        pub_folder = root.split('/')[-1]
        publication_keys = {}
        try:
            pub_data = json.load(open(root + '/publication.txt', 'r'))
            if 'url' in pub_data.keys():
                del pub_data['url']
            self.title = pub_data['title']
            self.authors = pub_data['authors']
            self.year = pub_data['year']
            if 'doi' not in pub_data:
                pub_data.update({'doi': None})
                print('ERROR: No doi')
                self.doi = None
            else:
                self.doi = pub_data['doi']
            if 'tags' not in pub_data:
                pub_data.update({'tags': None})
                print('ERROR: No tags')
                self.tags = None

            for key, value in pub_data.items():
                if isinstance(value, list):
                    value = json.dumps(value)
                else:
                    try:
                        value = int(value)
                    except BaseException:
                        pass

        except Exception as e:
            print(
                'ERROR: insufficient publication info {e}'.format(
                    **locals()))
            self.doi = None
            pub_data = {'title': None,
                        'authors': None,
                        'journal': None,
                        'volume': None,
                        'number': None,
                        'pages': None,
                        'year': None,
                        'publisher': None,
                        'doi': None,
                        'tags': None
                        }

        try:
            self.energy_corrections = json.load(
                open(root + '/energy_corrections.txt', 'r'))
        except BaseException:
            self.energy_corrections = {}

        if pub_data['title'] is None:
            self.title = root.split('/')[-1]
            pub_data.update({'title': self.title})
        if pub_data['authors'] is None:
            self.authors = [self.user]
            pub_data.update({'authors': self.authors})
        if pub_data['year'] is None:
            self.year = date.today().year
            pub_data.update({'year': self.year})

        self.pub_id = self.authors[0].split(',')[0].split(' ')[0] + \
            self.title.split(' ')[0].split('_')[0] + \
            str(self.year)

        self.cathub_db = '{}{}.db'.format(self.data_base, self.pub_id)

        pub_data.update({'pub_id': self.pub_id})
        self.pid = self.write_publication(pub_data)

    def read_gas(self, root):
        files = [f for f in os.listdir(root) if os.path.isfile(root + '/' + f)]
        traj_files = ['{}/{}'.format(root, f)
                      for f in files if f.endswith('.traj')]

        self.ase_ids_gas = {}
        self.traj_gas = {}

        for traj in traj_files:
            ase_id = None
            found = False

            if not ase_tools.check_traj(traj, self.strict, False):
                return

            chemical_composition = \
                ''.join(sorted(ase_tools.get_chemical_formula(
                    traj, mode='all')))
            chemical_composition_hill = ase_tools.get_chemical_formula(
                traj, mode='hill')
            energy = ase_tools.get_energies([traj])
            key_value_pairs = {"name": chemical_composition_hill,
                               'state': 'gas',
                               'epot': energy}

            id, ase_id = ase_tools.check_in_ase(
                traj, self.cathub_db)

            if ase_id is None:
                ase_id = ase_tools.write_ase(traj, self.cathub_db,
                                             self.user,
                                             **key_value_pairs)
            elif self.update:
                ase_tools.update_ase(self.cathub_db,
                                     id, **key_value_pairs)

            self.ase_ids_gas.update({chemical_composition: ase_id})
            self.traj_gas.update({chemical_composition: traj})

    def read_bulk(self, root, files):

        self.metal, self.crystal = root.split('/')[-1].split('_', 1)

        print('------------------------------------------------------')
        print('                    Surface:  {}'.format(self.metal))
        print('------------------------------------------------------')

        self.ase_ids = {}
        traj_bulk = ['{}/{}'.format(root, f)
                     for f in files if f.endswith('.traj') and 'bulk' in f][0]

        ase_id = None

        if not ase_tools.check_traj(traj_bulk, self.strict, False):
            return

        energy = ase_tools.get_energies([traj_bulk])

        key_value_pairs = {"name": self.metal,
                           'state': 'bulk',
                           'epot': energy}

        id, ase_id = ase_tools.check_in_ase(
            traj_bulk, self.cathub_db)  # self.ase_db)
        if ase_id is None:
            ase_id = ase_tools.write_ase(traj_bulk, self.cathub_db,
                                         self.user, **key_value_pairs)
        elif self.update:
            ase_tools.update_ase(self.cathub_db, id, **key_value_pairs)

            self.ase_ids.update({'bulk' + self.crystal: ase_id})

    def read_slab(self, root, files):
        self.facet = root.split('/')[-1]
        self.ase_facet = 'x'.join(list(self.facet))

        self.empty_traj = [
            '{}/{}'.format(root, f) for f in files if f.endswith('.traj')
            and 'empty' in f][0]

        ase_id = None
        if not ase_tools.check_traj(self.empty_traj, self.strict, False):
            return

        energy = ase_tools.get_energies([self.empty_traj])
        key_value_pairs = {"name": self.metal,
                           'state': 'star',
                           'epot': energy}

        key_value_pairs.update({'species': ''})

        id, ase_id = ase_tools.check_in_ase(
            self.empty_traj, self.cathub_db)  # self.ase_db)

        if ase_id is None:
            ase_id = ase_tools.write_ase(self.empty_traj, self.cathub_db,
                                         self.user, **key_value_pairs)
        elif self.update:
            ase_tools.update_ase(self.cathub_db, id, **key_value_pairs)
        self.ase_ids.update({'star': ase_id})

    def read_reaction(self, root, files):
        folder_name = root.split('/')[-1]

        self.reaction, self.sites = ase_tools.get_reaction_from_folder(
            folder_name)  # reaction dict

        print('----------- REACTION:  {} --> {} --------------'
              .format('+'.join(self.reaction['reactants']),
                      '+'.join(self.reaction['products'])))

        self.reaction_atoms, self.prefactors, self.prefactors_TS, \
            self.states = ase_tools.get_reaction_atoms(self.reaction)

        """Create empty dictionaries"""
        r_empty = ['' for n in range(len(self.reaction_atoms['reactants']))]
        p_empty = ['' for n in range(len(self.reaction_atoms['products']))]
        self.traj_files = {'reactants': r_empty[:],
                           'products': p_empty[:]}

        key_value_pairs = {}

        """ Match reaction gas species with their traj file """
        for key, mollist in self.reaction_atoms.items():
            for i, molecule in enumerate(mollist):
                if self.states[key][i] == 'gas':
                    assert molecule in self.ase_ids_gas.keys()
                    self.traj_files[key][i] = self.traj_gas[molecule]
                    species = ase_tools.clear_prefactor(
                        self.reaction[key][i])
                    key_value_pairs.update(
                        {'species': ase_tools.clear_state(species)})
                    self.ase_ids.update({species: self.ase_ids_gas[molecule]})

    def read_energies(self, root, files):

        self.key_value_pairs_reaction = None
        if 'TS' in self.traj_files:
            del self.traj_files['TS']
        if 'TSempty' in self.traj_files:
            del self.traj_files['TSempty']

        traj_slabs = [f for f in files if f.endswith('.traj') and
                      'gas' not in f]

        if not self.debug:
            assert len(traj_slabs) > 0, \
                'Need at least one file in {}!'.format(root)
        else:
            try:
                assert len(traj_slabs) > 0
            except BaseException:
                print('Need at least one file in {}!'.format(root))
                return

        n_atoms = np.array([])
        ts_i = None
        tsempty_i = None
        chemical_composition_slabs = []
        breakloop = False
        for i, f in enumerate(traj_slabs):
            if 'empty' in f and 'TS' in f:
                tsempty_i = i
            elif 'TS' in f:
                ts_i = i

            traj = '{}/{}'.format(root, f)

            if not ase_tools.check_traj(traj, self.strict, False):
                return

            chemical_composition_slabs = \
                np.append(chemical_composition_slabs,
                          ase_tools.get_chemical_formula(traj, mode='all'))
            n_atoms = np.append(n_atoms, ase_tools.get_number_of_atoms(traj))

        traj_empty = self.empty_traj
        empty_atn = ase_tools.get_atomic_numbers(traj_empty)

        prefactor_scale = copy.deepcopy(self.prefactors)
        for key1, values in prefactor_scale.items():
            prefactor_scale[key1] = [1 for v in values]

        key_value_pairs = {}

        key_value_pairs.update({'name':
                                ase_tools.get_chemical_formula(traj_empty),
                                # 'site': self.sites,
                                'facet': self.ase_facet,
                                'layers': ase_tools.get_n_layers(traj_empty),
                                'state': 'star'})

        """ Write empty slab to ASE"""
        ase_id = None
        id, ase_id = ase_tools.check_in_ase(traj_empty, self.cathub_db)
        for key, mollist in self.reaction_atoms.items():
            if '' in mollist:
                n = mollist.index('')
                self.traj_files[key][n] = traj_empty
                key_value_pairs.update({'species': ''})
            if ase_id is None:
                ase_id = ase_tools.write_ase(traj_empty, self.cathub_db,
                                             self.user, **key_value_pairs)
            elif self.update:
                ase_tools.update_ase(self.cathub_db, id, **key_value_pairs)

            self.ase_ids.update({'star': ase_id})

        """ Handle other slabs"""
        for i, f in enumerate(traj_slabs):
            traj = '{}/{}'.format(root, f)
            atns = ase_tools.get_atomic_numbers(traj)
            if not (np.array(atns) > 8).any() and \
               (np.array(empty_atn) > 8).any():
                print("Only molecular species in traj file: {}".format(traj))
                continue

            # Get supercell size relative to empty slab
            supercell_factor = 1
            if len(atns) > len(empty_atn) * 2:  # different supercells
                supercell_factor = len(res_slab_atn) // len(empty_atn)

            # Atomic numbers of adsorbate
            ads_atn = atns.copy()
            for atn in empty_atn * supercell_factor:
                ads_atn.remove(atn)
            ads_atn = sorted(ads_atn)

            ase_id = None
            id, ase_id = ase_tools.check_in_ase(traj, self.cathub_db)
            key_value_pairs.update({'epot': ase_tools.get_energies([traj])})

            if i == ts_i:  # transition state
                self.traj_files.update({'TS': [traj]})
                self.prefactors.update({'TS': [1]})
                prefactor_scale.update({'TS': [1]})
                key_value_pairs.update({'species': 'TS'})
                if ase_id is None:
                    ase_id = ase_tools.write_ase(traj, self.cathub_db,
                                                 self.user, **key_value_pairs)
                elif self.update:
                    ase_tools.update_ase(self.cathub_db, id, **key_value_pairs)
                self.ase_ids.update({'TSstar': ase_id})
                continue

            if i == tsempty_i:  # empty slab for transition state
                self.traj_files.update({'TSempty': [traj]})
                self.prefactors.update({'TSempty': [1]})
                prefactor_scale.update({'TSempty': [1]})
                key_value_pairs.update({'species': ''})
                if ase_id is None:
                    ase_id = ase_tools.write_ase(traj, self.cathub_db,
                                                 self.user, **key_value_pairs)
                elif self.update:
                    ase_tools.update_ase(self.cathub_db, id, **key_value_pairs)
                self.ase_ids.update({'TSemptystar': ase_id})
                continue

            found = False
            for key, mollist in self.reaction_atoms.items():
                if found:
                    break
                for n, molecule in enumerate(mollist):
                    if found:
                        break
                    molecule_atn = ase_tools.get_numbers_from_formula(molecule)
                    for n_ads in range(1, 5):
                        mol_atn = sorted(molecule_atn * n_ads)
                        if ads_atn == mol_atn and \
                           self.states[key][n] == 'star':
                            self.traj_files[key][n] = traj
                            species = ase_tools.clear_prefactor(
                                self.reaction[key][n])
                            id, ase_id = ase_tools.check_in_ase(
                                traj, self.cathub_db)
                            key_value_pairs.update(
                                {'species':
                                 ase_tools.clear_state(
                                     species),
                                 'n': n_ads,
                                 'site': self.sites[species]})
                            if ase_id is None:
                                ase_id = ase_tools.write_ase(
                                    traj, self.cathub_db, self.user,
                                    **key_value_pairs)
                            elif self.update:
                                ase_tools.update_ase(
                                    self.cathub_db, id, **key_value_pairs)
                            self.ase_ids.update({species: ase_id})
                            found = True
                            break

            if n_ads > 1:
                for key1, values in prefactor_scale.items():
                    for mol_i in range(len(values)):
                        if self.states[key1][mol_i] == 'gas':
                            prefactor_scale[key1][mol_i] = n_ads

            if supercell_factor > 1:
                for key2, values in prefactor_scale.items():
                    for mol_i in range(len(values)):
                        if self.reaction[key2][mol_i] == 'star':
                            prefactor_scale[key2][mol_i] *= supercell_factor

        surface_composition = self.metal
        chemical_composition = ase_tools.get_chemical_formula(traj_empty)

        prefactors_final = copy.deepcopy(self.prefactors)
        for key in self.prefactors:
            for i, v in enumerate(self.prefactors[key]):
                prefactors_final[key][i] = self.prefactors[key][i] * \
                    prefactor_scale[key][i]

        reaction_energy = None
        activation_energy = None
        try:
            reaction_energy, activation_energy = \
                ase_tools.get_reaction_energy(
                    self.traj_files, self.reaction,
                    self.reaction_atoms,
                    self.states, prefactors_final,
                    self.prefactors_TS,
                    self.energy_corrections)

        except BaseException:
            if self.debug:
                print('ERROR: reaction energy failed for files in: {}'
                      .format(root))
            else:
                raise RuntimeError(
                    'Reaction energy failed for files in: {}'.format(root))

        expr = -10 < reaction_energy < 10

        if not ase_tools.debug_assert(
                expr, 'reaction energy is wrong: {} eV: {}'
                .format(reaction_energy, root),
                self.debug):
            return

        expr = activation_energy is None \
               or reaction_energy < activation_energy < 5
        if not ase_tools.debug_assert(expr,
                                      'activation energy is wrong: {} eV: {}'
                                      .format(activation_energy, root),
                                      self.debug):
            print(self.traj_files, prefactors_final, self.prefactors_TS)

        reaction_info = {'reactants': {},
                         'products': {}}

        for key in ['reactants', 'products']:
            for i, r in enumerate(self.reaction[key]):
                r = ase_tools.clear_prefactor(r)
                reaction_info[key].update({r: self.prefactors[key][i]})

        self.key_value_pairs_reaction = {
            'chemical_composition': chemical_composition,
            'surface_composition': surface_composition,
            'facet': self.facet,
            'sites': self.sites,
            'coverages': self.coverages,
            'reactants': reaction_info['reactants'],
            'products': reaction_info['products'],
            'reaction_energy': reaction_energy,
            'activation_energy': activation_energy,
            'dft_code': self.DFT_code,
            'dft_functional': self.DFT_functional,
            'pub_id': self.pub_id,
            'doi': self.doi,
            'year': int(self.year),
            'ase_ids': self.ase_ids,
            'energy_corrections': self.energy_corrections,
            'username': self.user}


def read_name_from_folder(root):
    folder_name = root.split('/')[-1]
    return folder_name
