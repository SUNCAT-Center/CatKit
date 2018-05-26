from ase.db.sqlite import SQLite3Database
import sqlite3
import json


init_commands = [
    """ CREATE TABLE publication (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pub_id text UNIQUE,
    title text,
    authors text,
    journal text,
    volume text,
    number text,
    pages text,
    year integer,
    publisher text,
    doi text,
    tags text
    );""",

    """CREATE TABLE publication_system (
    ase_id text REFERENCES systems(unique_id),
    pub_id text REFERENCES publication(pub_id),
    PRIMARY KEY (pub_id, ase_id)
    );""",

    """CREATE TABLE reaction (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    chemical_composition text,
    surface_composition text,
    facet text,
    sites text,
    coverages text,
    reactants text,
    products text,
    reaction_energy real,
    activation_energy real,
    dft_code text,
    dft_functional text,
    username text,
    pub_id text,
    FOREIGN KEY (pub_id) REFERENCES publication(pub_id)
    );""",

    """ CREATE TABLE reaction_system (
    name text,
    energy_correction real,
    ase_id text,
    id integer,
    FOREIGN KEY (ase_id) REFERENCES systems(unique_id),
    FOREIGN KEY (id) REFERENCES reaction(id)
    );"""]


class CathubSQLite:
    """TODO: Please document me."""

    def __init__(self, filename):
        """TODO: Please document me.

        Parameters
        ----------
        default : str
            used for autoincrement id
        """
        self.filename = filename
        self.initialized = False
        self.default = 'NULL'
        self.connection = None
        self.id = None
        self.pid = None

    def _connect(self):
        return sqlite3.connect(self.filename, timeout=600)

    def __enter__(self):
        assert self.connection is None
        self.connection = self._connect()
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if exc_type is None:
            self.connection.commit()
        else:
            self.connection.rollback()
        self.connection.close()
        self.connection = None

    def _initialize(self, con):
        if self.initialized:
            return
        SQLite3Database()._initialize(con)

        self._metadata = {}

        cur = con.execute(
            'SELECT COUNT(*) FROM sqlite_master WHERE name="reaction"')

        if cur.fetchone()[0] == 0:
            for init_command in init_commands:
                con.execute(init_command)
            self.id = 1
            self.pid = 1
            con.commit()

        if self.id is None:
            self.id = self.get_last_id(cur) + 1
            self.pid = self.get_last_pub_id(cur) + 1

        self.initialized = True

    def read(self, id, table='reaction'):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        cur.execute('SELECT * FROM \n {} \n WHERE \n {}.id={}'.format(
            table, table, id))
        row = cur.fetchall()

        if len(row) == 14:  # Old schema
            row = row.insert(5, 'None')

        return row

    def write_publication(self, values):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        pid = self.pid

        values = (pid,
                  values['pub_id'],
                  values['title'],
                  json.dumps(values['authors']),
                  values['journal'],
                  values['volume'],
                  values['number'],
                  values['pages'],
                  values['year'],
                  values['publisher'],
                  values['doi'],
                  json.dumps(values['tags']))

        q = ', '.join('?' * len(values))
        cur.execute('INSERT OR IGNORE INTO publication VALUES ({})'.format(q),
                    values)
        self.pid = pid + 1

        return pid

    def check_ase_ids(self, values, ase_ids):
        ase_values = ase_ids.values()
        assert len(set(ase_values)) == len(ase_values), 'Duplicate ASE ids!'

        reaction_species = set(list(values['reactants'].keys()) +
                               list(values['products'].keys()))

        n_split = 0
        for spec in ase_ids.keys():
            if '_' in spec:
                n_split += 1

        assert len(reaction_species) <= len(ase_values) + n_split, \
            'ASE ids missing!'
        return

    def write(self, values, data=None):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        id = self.id

        pub_id = values['pub_id']
        ase_ids = values['ase_ids']
        energy_corrections = values['energy_corrections']

        if ase_ids is not None:
            self.check_ase_ids(values, ase_ids)
        else:
            ase_ids = {}
        values = (id,
                  values['chemical_composition'],
                  values['surface_composition'],
                  values['facet'],
                  json.dumps(values['sites']),
                  json.dumps(values['coverages']),
                  json.dumps(values['reactants']),
                  json.dumps(values['products']),
                  values['reaction_energy'],
                  values['activation_energy'],
                  values['dft_code'],
                  values['dft_functional'],
                  values['username'],
                  values['pub_id']
                  )

        q = ', '.join('?' * len(values))
        cur.execute('INSERT INTO reaction VALUES ({})'.format(q),
                    values)
        reaction_structure_values = []

        for name, ase_id in ase_ids.items():
            if name in energy_corrections:
                energy_correction = energy_corrections[name]
            else:
                energy_correction = 0
            reaction_structure_values.append([name, energy_correction,
                                              ase_id, id])
            insert_statement = """INSERT OR IGNORE INTO
            publication_system(ase_id, pub_id) VALUES (?, ?)"""
            cur.execute(insert_statement, [ase_id, pub_id])

        cur.executemany('INSERT INTO reaction_system VALUES (?, ?, ?, ?)',
                        reaction_structure_values)

        if self.connection is None:
            con.commit()
            con.close()

        self.id = id + 1

        return id

    def update(self, id, values, key_names='all'):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        pub_id = values['pub_id']
        ase_ids = values['ase_ids']
        energy_corrections = values['energy_corrections']

        if ase_ids is not None:
            self.check_ase_ids(values, ase_ids)
        else:
            ase_ids = {}

        key_list, value_list = get_key_value_list(key_names, values)
        N_keys = len(key_list)

        value_strlist = get_value_strlist(value_list)
        execute_str = ', '.join('{}={}'.format(key_list[i], value_strlist[i])
                                for i in range(N_keys))

        update_command = 'UPDATE reaction SET {} WHERE id = {};'\
            .format(execute_str, id)

        cur.execute(update_command)

        delete_command = 'DELETE from reaction_system WHERE id = {}'.format(id)
        cur.execute(delete_command)

        reaction_structure_values = []
        for name, ase_id in ase_ids.items():
            reaction_structure_values.append([name,
                                              energy_corrections.get(name),
                                              ase_id, id])
            insert_statement = """INSERT OR IGNORE INTO
            publication_system(ase_id, pub_id) VALUES (?, ?)"""
            cur.execute(insert_statement, [ase_id, pub_id])

        cur.executemany('INSERT INTO reaction_system VALUES (?, ?, ?, ?)',
                        reaction_structure_values)

        if self.connection is None:
            con.commit()
            con.close()
        return id

    def get_last_id(self, cur):
        cur.execute('SELECT seq FROM sqlite_sequence WHERE name="reaction"')
        result = cur.fetchone()
        if result is not None:
            id = result[0]
        else:
            id = 0
        return id

    def get_last_pub_id(self, cur):
        cur.execute('SELECT seq FROM sqlite_sequence WHERE name="publication"')
        result = cur.fetchone()
        if result is not None:
            pid = result[0]
        else:
            pid = 0
        return pid

    def check(self, chemical_composition, reaction_energy):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        statement = """SELECT reaction.id FROM reaction WHERE
        reaction.chemical_composition=? and reaction.reaction_energy=?"""
        argument = [chemical_composition, reaction_energy]

        cur.execute(statement, argument)
        rows = cur.fetchall()
        if len(rows) > 0:
            id = rows[0][0]
        else:
            id = None
        return id

    def check_reaction_on_surface(self, chemical_composition, reactants,
                                  products):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        statement = """SELECT reaction.id FROM reaction WHERE
        reaction.chemical_composition='{}' and reaction.reactants='{}'
        and reaction.products='{}';""".format(chemical_composition,
                                              json.dumps(reactants),
                                              json.dumps(products))

        cur.execute(statement)
        rows = cur.fetchall()
        if len(rows) > 0:
            id = rows[0][0]
        else:
            id = None
        return id

    def check_publication(self, pub_id):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        statement = """
        SELECT id FROM publication WHERE publication.pub_id=?"""
        argument = [pub_id]

        cur.execute(statement, argument)
        rows = cur.fetchall()
        if len(rows) > 0:
            id = rows[0][0]
        else:
            id = None
        return id

    def check_publication_structure(self, pub_id, ase_id):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        statement = """
        SELECT id FROM publication_system WHERE
        publication.pub_id=? and publication.ase_id=?"""
        argument = [pub_id, ase_id]

        cur.execute(statement, argument)
        rows = cur.fetchall()
        if len(rows) > 0:
            id = rows[0][0]
        else:
            id = None
        return id


def get_key_value_str(values):
    key_str = """chemical_composition, surface_composition, facet, sites,
    reactants,products, reaction_energy, activation_energy, dft_code,
    dft_functional, publication, doi, year, ase_ids, user"""
    value_str = "'{}'".format(values[1])
    for v in values[2:]:
        try:
            # Unicode in python2 - must be a better way of handling this.
            if isinstance(v, unicode):
                v = v.encode('ascii', 'ignore')
        except NameError:
            pass
        if isinstance(v, str) or isinstance(v, dict):
            value_str += ", '{}'".format(v)
        elif v is None or v == '':
            value_str += ", {}".format('NULL')
        else:
            value_str += ", {}".format(v)

    return key_str, value_str


def get_key_value_list(key_list, values, table='reaction'):
    total_keys = {'reaction': ['chemical_composition', 'surface_composition',
                               'facet', 'sites', 'coverages', 'reactants',
                               'products', 'reaction_energy',
                               'activation_energy', 'dft_code',
                               'dft_functional', 'username', 'pub_id'],
                  'publication': ['pub_id', 'title', 'authors', 'journal',
                                  'volume', 'number', 'pages', 'year',
                                  'publisher', 'doi', 'tags'],
                  'reaction_system': ['name', 'energy_correction',
                                      'ase_id', 'reaction_id'],
                  'publication_system': ['ase_id, pub_id']}
    total_key_list = total_keys[table]

    if key_list == 'all':
        key_list = total_key_list
    else:
        for key in key_list:
            assert key in total_key_list

    value_list = [values[key] for key in key_list]
    return key_list, value_list


def get_value_strlist(value_list):
    value_strlist = []
    for v in value_list:
        try:
            # Unicode in python2 - must be a better way of handling this.
            if isinstance(v, unicode):
                v = v.encode('ascii', 'ignore')
        except NameError:
            pass
        if isinstance(v, dict):
            v = json.dumps(v)
            value_strlist.append("'{}'".format(v))
        elif isinstance(v, str):
            value_strlist.append("'{}'".format(v))
        elif v is None or v == '':
            value_strlist.append("{}".format('NULL'))
        else:
            value_strlist.append("{}".format(v))

    return value_strlist
