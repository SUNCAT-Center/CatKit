from .cathubsqlite import CathubSQLite
from .tools import get_bases
from . import ase_tools
import numpy as np
import os
import copy
import json


class FolderReader:
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
        self.reaction_level = 4
        self.metal_level = 5
        self.facet_level = 6
        self.site_level = None
        self.final_level = 6
        self.title = None
        self.authors = None
        self.omit_folders = []

        self.coverages = None
        # print(self.site_level)
        # if os.path.isfile(user_file):
        #    user_spec = json.load(open(user_file, 'r'))
        #    self.__dict__.update(user_spec)

    def read(self, skip=[], goto_reaction=None):
        if len(skip) > 0:
            for skip_f in skip:
                self.omit_folders.append(skip_f)
        up = 0
        found_reaction = False
        for root, dirs, files in os.walk(self.user_base):
            for omit_folder in self.omit_folders:  # user specified omit_folder
                if omit_folder in dirs:
                    dirs.remove(omit_folder)
            level = len(root.split("/")) - self.user_base_level

            if level == self.pub_level:
                self.read_pub(root)

            if level == self.DFT_level:
                self.DFT_code = self.read_name_from_folder(root)

            if level == self.XC_level:
                self.DFT_functional = self.read_name_from_folder(root)

            if level == self.reaction_level:
                # print(root.split("/")[-1])
                if goto_reaction is not None and not found_reaction:
                    if not root.split("/")[-1] == goto_reaction:
                        dirs[:] = []
                        continue
                    else:
                        found_reaction = True

                self.read_reaction(root, files)

            if level == self.metal_level:
                up = 0

            if level == self.metal_level + up:
                self.read_metal(root)
                # if self.user == 'roling':
                #    if self.metal ==\
                #       self.reaction['reactants'][0].replace('star', ''):
                #        up = 1
                #        continue

            if level == self.facet_level + up \
               and not level == self.metal_level + up:
                self.read_facet(root)

            if self.site_level is not None and self.site_level != 'None':
                if level == int(self.site_level) + up\
                   and not level == self.facet_level + up:
                    self.read_site(root)
                    try:

                        self.sites = int(self.sites)
                    except BaseException:
                        pass

            else:
                self.sites = ''

            if level == self.final_level + up:
                self.read_final(root, files)
                if self.key_value_pairs_reaction is not None:
                    yield self.key_value_pairs_reaction

    def write(self, skip=[], goto_reaction=None):
        for key_values in self.read(skip=skip, goto_reaction=goto_reaction):
            with CathubSQLite(self.cathub_db) as db:
                id = db.check(
                    key_values['chemical_composition'],
                    key_values['reaction_energy'])
                #print('Allready in reaction db with row id = {}'.format(id))
                if id is None:
                    id = db.write(key_values)
                    print('Written to reaction db row id = {}'.format(id))
                elif self.update:
                    db.update(id, key_values)
                    print('Updated reaction db row id = {}'.format(id))
                else:
                    print('Allready in reaction db with row id = {}'.format(id))

    def write_publication(self, pub_data):
        with CathubSQLite(self.cathub_db) as db:
            pid = db.check_publication(self.pub_id)
            if pid is None:
                pid = db.write_publication(pub_data)
                print('Written to publications db row id = {}'.format(pid))
            else:
                print('Allready in reaction db with row id = {}'.format(pid))
        return pid

    def update_sqlite(self, skip=[], goto_reaction=None, key_names='all'):
        for key_values in self.read(skip=skip, goto_reaction=goto_reaction):
            with CathubSQLite(self.cathub_db) as db:
                id = db.check(key_values['reaction_energy'])
                #print('Allready in reaction db with row id = {}'.format(id))
                if id is not None:
                    db.update(id, key_values, key_names)

    def read_name_from_folder(self, root):
        folder_name = root.split('/')[-1]
        return folder_name

    def read_pub(self, root):
        pub_folder = root.split('/')[-1]
        #self.ase_db = '{}atoms_{}.db'.format(self.data_base, pub_folder)
        #self.catapp_db = '{}catapp_{}.db'.format(self.data_base, pub_folder)
        # assert 'publication.txt' in files
        publication_keys = {}
        try:
            pub_data = json.load(open(root + '/publication.txt', 'r'))
            if 'url' in pub_data.keys():
                del pub_data['url']
            self.reference = json.dumps(pub_data)
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
               # publication_keys.update({'publication_' + key: value})

        except BaseException:
            print('ERROR: insufficient publication info')
            self.doi = None
            pub_data = {'title': self.title,
                        'authors': self.authors,
                        'journal': None,
                        'volume': None,
                        'number': None,
                        'pages': None,
                        'year': None,
                        'publisher': None,
                        'doi': None,
                        'tags': None
                        }
            self.reference = json.dumps(pub_data)

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
            self.year = 2018
            pub_data.update({'year': self.year})

        self.pub_id = self.authors[0].split(',')[0].split(' ')[0] + \
            self.title.split(' ')[0].split('_')[0] + \
            str(self.year)

        self.cathub_db = '{}{}.db'.format(self.data_base, self.pub_id)

        pub_data.update({'pub_id': self.pub_id})
        self.pid = self.write_publication(pub_data)

        #self.publication_keys = publication_keys

    def read_reaction(self, root, files):
        folder_name = root.split('/')[-1]

        # try:
        self.reaction, sites = ase_tools.get_reaction_from_folder(
            folder_name)  # reaction dict
        # except:
        #    print('ERROR: omitting directory {}'.format(root))
        #    dirs = []
        #    return False

        print('-------------- REACTION:  {} --> {} -----------------'
              .format('+'.join(self.reaction['reactants']),
                      '+'.join(self.reaction['products'])))

        self.reaction_atoms, self.prefactors, self.prefactors_TS, \
            self.states = ase_tools.get_reaction_atoms(self.reaction)

        # Create empty dictionaries
        r_empty = ['' for n in range(len(self.reaction_atoms['reactants']))]
        p_empty = ['' for n in range(len(self.reaction_atoms['products']))]
        self.traj_files = {'reactants': r_empty[:],
                           'products': p_empty[:]}

        self.ase_ids = {}

        traj_gas = ['{}/{}'.format(root, f)
                    for f in files if f.endswith('.traj')]

        if len(traj_gas) == 0 and \
                np.any(['gas' in self.states.values()[s]
                        for s in range(len(self.states.values()))]):

            traj_gas = []
            for root2, dirs2, files2 in os.walk(root):
                for f2 in [f2 for f2 in files2 if f2.endswith('gas.traj')]:
                    traj_gas.append('{}/{}'.format(root2, f2))

        key_value_pairs = {}  # copy.deepcopy(self.publication_keys)
        for traj in traj_gas:
            ase_id = None
            found = False

            if not ase_tools.check_traj(traj, self.strict, False):
                return
            chemical_composition = ''.join(
                sorted(
                    ase_tools.get_chemical_formula(
                        traj, mode='all')))
            chemical_composition_hill = ase_tools.get_chemical_formula(
                traj, mode='hill')

            energy = ase_tools.get_energies([traj])
            key_value_pairs.update({"name": chemical_composition_hill,
                                    'state': 'gas',
                                    'epot': energy})

            id, ase_id = ase_tools.check_in_ase(
                traj, self.cathub_db)  # self.ase_db)

            for key, mollist in self.reaction_atoms.items():
                for i, molecule in enumerate(mollist):
                    if molecule == chemical_composition \
                            and self.states[key][i] == 'gas':
                        # Should only be found once?
                        assert found is False, \
                            root + ' ' + chemical_composition
                        found = True
                        self.traj_files[key][i] = traj

                        species = ase_tools.clear_prefactor(
                            self.reaction[key][i])
                        key_value_pairs.update(
                            {'species': ase_tools.clear_state(species)})
                        if ase_id is None:
                            ase_id = ase_tools.write_ase(
                                traj, self.cathub_db,
                                self.user, **key_value_pairs)
                        elif self.update:
                            ase_tools.update_ase(
                                self.cathub_db, id, **key_value_pairs)
                        self.ase_ids.update({species: ase_id})

            # if found is False:
            #    print('{} file is not part of reaction, include as reference'\
            #        .format(f))
            #    self.ase_ids.update({chemical_composition_hill + 'gas': ase_id})

    def read_metal(self, root):
        self.metal = root.split('/')[-1]
        if self.metal_level == self.facet_level:
            if len(self.metal.split('_')) == 2:
                self.metal, self.facet = self.metal.split('_')
                self.sites = ''
                self.ase_facet = self.facet
                if 'x' in self.facet and not '-' not in self.facet:
                    facetstr = 'x'.join('{}' for f in self.facet)
                    self.ase_facet = facetstr.format(*self.facet)

            else:
                self.facet = None
                self.ase_facet = ''
                self.sites = ''
        print('--------------- METAL: {} ---------------'.format(self.metal))

    def read_facet(self, root):
        folder_name = root.split('/')[-1]

        if not self.facet_level == self.site_level or self.site_level is None:
            self.facet = folder_name
        else:
            split = folder_name.split('_')
            if len(split) == 1:
                split = split[0].split('-')
            self.facet, self.sites = split
        self.ase_facet = self.facet
        if 'x' in self.facet and not '-' not in self.facet:
            facetstr = 'x'.join('{}' for f in self.facet)
            self.ase_facet = facetstr.format(*self.facet)
        print('--------------- FACET: {} ---------------'.format(self.facet))

    def read_site(self, root):
        dirjoin = '_'.join(info for info in root.split('/')
                           [self.site_level + self.user_base_level - 1:])
        self.sites = dirjoin

    def read_final(self, root, files):
        self.key_value_pairs_reaction = None
        if 'TS' in self.traj_files:
            del self.traj_files['TS']
        if 'TSempty' in self.traj_files:
            del self.traj_files['TSempty']
        traj_slabs = [f for f in files if f.endswith('.traj') and
                      'gas' not in f]
        # if traj_slabs == []:
        #    return
        if not self.debug:
            assert len(traj_slabs) > 1, \
                'Need at least two files in {}!'.format(root)
        else:
            try:
                assert len(traj_slabs) > 1, \
                    'Need at least two files in {}!'.format(root)
            except BaseException:
                print('Need at least two files in {}!'.format(root))
                return

        n_atoms = np.array([])
        empty_i = None
        ts_i = None
        tsempty_i = None
        chemical_composition_slabs = []
        breakloop = False
        for i, f in enumerate(traj_slabs):
            if 'empty' in f and 'TS' in f:
                tsempty_i = i
            elif 'empty' in f:
                empty_i = i
            elif 'TS' in f:
                ts_i = i

            traj = '{}/{}'.format(root, f)
            if not ase_tools.check_traj(traj, self.strict, False):
                breakloop = True
                break
            chemical_composition_slabs = \
                np.append(chemical_composition_slabs,
                          ase_tools.get_chemical_formula(traj, mode='all'))
            n_atoms = np.append(n_atoms, ase_tools.get_number_of_atoms(traj))

        if breakloop:
            return

        # Empty slab has least atoms
        if empty_i is None:
            empty_i = np.argmin(n_atoms)
        traj_empty = root + '/' + traj_slabs[empty_i]

        empty_atn = ase_tools.get_atomic_numbers(traj_empty)

        # Identify TS
        # if ts_i is not None:
        #    traj_TS = root + '/' + traj_slabs[ts_i]
        #    self.traj_files.update({'TS': [traj_TS]})
        #    self.prefactors.update({'TS': [1]})
        #    TS_id = {get_chemical_formula(traj_TS): ase_id}

        # elif ts_i is None and len(traj_slabs) > len(reaction) + 1:
        #raise AssertionError('which one is the transition state???')
        # else:
        #    TS_id = None
        #    activation_energy = None

        prefactor_scale = copy.deepcopy(self.prefactors)
        for key1, values in prefactor_scale.items():
            prefactor_scale[key1] = [1 for v in values]

        #prefactor_scale_ads = copy.deepcopy(prefactor_scale)

        key_value_pairs = {}  # copy.deepcopy(self.publication_keys)
        key_value_pairs.update(
            {
                'name': ase_tools.get_chemical_formula(traj_empty),
                'site': self.sites,
                'facet': self.ase_facet,
                'layers': ase_tools.get_n_layers(traj_empty),
                'state': 'star'})

        for i, f in enumerate(traj_slabs):
            traj = '{}/{}'.format(root, f)
            ase_id = None
            id, ase_id = ase_tools.check_in_ase(traj, self.cathub_db)
            found = False
            key_value_pairs.update({'epot': ase_tools.get_energies([traj])})

            if i == ts_i:
                found = True
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

            if i == tsempty_i:
                found = True
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

            elif i == empty_i:
                found = True
                for key, mollist in self.reaction_atoms.items():
                    if '' in mollist:
                        n = mollist.index('')
                        self.traj_files[key][n] = traj
                key_value_pairs.update({'species': ''})
                if ase_id is None:
                    ase_id = ase_tools.write_ase(traj, self.cathub_db,
                                                 self.user, **key_value_pairs)
                elif self.update:
                    ase_tools.update_ase(self.cathub_db, id, **key_value_pairs)
                self.ase_ids.update({'star': ase_id})
                continue

            res_atn = ase_tools.get_atomic_numbers(traj)
            if not (
                    np.array(res_atn) > 8).any() and (
                    np.array(empty_atn) > 8).any():
                print("Only molecular species in traj file: {}".format(traj))
                continue

            for atn in empty_atn:
                res_atn.remove(atn)

            res_atn = sorted(res_atn)  # residual atomic numbers

            supercell_factor = 1
            n_ads = 1
            for key, mollist in self.reaction_atoms.items():
                if found:
                    continue
                for n, molecule in enumerate(mollist):
                    if found:
                        continue
                    molecule_atn = ase_tools.get_numbers_from_formula(molecule)
                    for k in range(1, 5):
                        if found:
                            continue
                        mol_atn = sorted(molecule_atn * k)
                        if res_atn == mol_atn and self.states[key][n] == 'star':
                            found = True
                        elif len(res_atn) >= len(empty_atn):
                            res_slab_atn = res_atn[:]
                            for atn in mol_atn:
                                if atn in res_slab_atn:
                                    res_slab_atn.remove(atn)
                            try:
                                supercell_factor = int(
                                    len(res_slab_atn) / len(empty_atn))
                                if sorted(empty_atn * supercell_factor) \
                                        == sorted(res_slab_atn):
                                    found = True
                            except BaseException:
                                continue

                        if found:
                            n_ads = k
                            self.traj_files[key][n] = traj

                            species = ase_tools.clear_prefactor(
                                self.reaction[key][n])

                            id, ase_id = ase_tools.check_in_ase(
                                traj, self.cathub_db)
                            key_value_pairs.update({'species':
                                                    ase_tools.clear_state(species),
                                                    'n': n_ads})
                            if ase_id is None:
                                ase_id = ase_tools.write_ase(
                                    traj, self.cathub_db, self.user, **key_value_pairs)
                            elif self.update:
                                ase_tools.update_ase(
                                    self.cathub_db, id, **key_value_pairs)
                            self.ase_ids.update({species: ase_id})

            if n_ads > 1:
                for key1, values in prefactor_scale.items():
                    for mol_i in range(len(values)):
                        #prefactor_scale_ads[key1][mol_i] = n_ads
                        if self.states[key1][mol_i] == 'gas':
                            prefactor_scale[key1][mol_i] = n_ads

            if supercell_factor > 1:
                for key2, values in prefactor_scale.items():
                    for mol_i in range(len(values)):
                        if self.reaction[key2][mol_i] == 'star':
                            prefactor_scale[key2][mol_i] *= supercell_factor + 1

        # Transition state has higher energy
        # if len(np.unique(chemical_compositions)) > len(chemical_compositions):
        #    for chemical_composition in chemical_compositions:

        surface_composition = self.metal  # get_surface_composition(traj_empty)
        chemical_composition = ase_tools.get_chemical_formula(traj_empty)

        prefactors_final = copy.deepcopy(self.prefactors)
        for key in self.prefactors:
            for i, v in enumerate(self.prefactors[key]):
                prefactors_final[key][i] = self.prefactors[key][i] * \
                    prefactor_scale[key][i]

        reaction_energy = None
        activation_energy = None

        if self.user == 'jschum':
            if 'H2gas' in self.ase_ids.keys():
                self.energy_corrections.update({'H2gas': 0.1})
            if 'CH3CHOgas' in self.ase_ids.keys():
                self.energy_corrections.update({'CH3CHOgas': 0.15})

            if self.ase_ids.get('TSemptystar') == self.ase_ids.get('star'):
                del self.ase_ids['TSemptystar']

        try:
            reaction_energy, activation_energy = \
                ase_tools.get_reaction_energy(
                    self.traj_files,
                    self.reaction,
                    self.reaction_atoms,
                    self.states,
                    prefactors_final,
                    self.prefactors_TS,
                    self.energy_corrections)

        except BaseException:
            if self.debug:
                print('ERROR: reaction energy failed for files in: {}'.format(root))
                #flag = True
            else:
                raise RuntimeError(
                    'Reaction energy failed for files in: {}'.format(root))

        expr = -10 < reaction_energy < 10

        if not ase_tools.debug_assert(expr,
                                      'reaction energy is wrong: {} eV: {}'
                                      .format(reaction_energy, root),
                                      self.debug):

            return

        expr = ase_tools.activation_energy is None or reaction_energy < activation_energy < 5
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
            'year': int(
                self.year),
            'ase_ids': self.ase_ids,
            'energy_corrections': self.energy_corrections,
            'username': self.user}
