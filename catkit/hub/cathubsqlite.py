import sys
from ase.db.sqlite import SQLite3Database
import sqlite3
import json
from past.utils import PY2
from tabulate import tabulate


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
    """Class for managing SQLite3 database for reaction energies,
    publications and atomic structures. Builds on top of the ASE database for
    atomic strucutres https://wiki.fysik.dtu.dk/ase/ase/db/db.html with four
    additional tables:

       publication: publication info

       publication_system: one-to-many mapping between publication table and
           systems table in ASE database

       reaction: reaction energies for surfaces

       reaction_system: mamy-to-many mapping between reaction table and
           systems table in ASE database

    Connect to a database object:

        db = CathubSQLite('yourdbfile.db')

    Set up a connection for several manipulations:

        with db as CathubSQLite('yourdbfile.db'):
            Do your work...

    Parameters
    ----------
    filename : str
        name of database file
    """

    def __init__(self, filename, stdin=sys.stdin, stdout=sys.stdout):

        assert filename.endswith('.db'), 'filename should have .db extension'
        self.filename = filename
        self.initialized = False
        self.default = 'NULL'
        self.connection = None
        self.stdin = stdin
        self.stdout = stdout

    def _connect(self):
        return sqlite3.connect(self.filename, timeout=600)

    def __enter__(self):
        """Set connection upon entry using with statement"""
        assert self.connection is None
        self.connection = self._connect()
        return self

    def __exit__(self, exc_type, exc_value, tb):
        """Commit changes upon exit"""
        if exc_type is None:
            self.connection.commit()
        else:
            self.connection.rollback()
        self.connection.close()
        self.connection = None

    def _initialize(self, con):
        """Set up tables in SQL"""
        if self.initialized:
            return

        SQLite3Database()._initialize(con)  # ASE db initialization

        cur = con.execute(
            'SELECT COUNT(*) FROM sqlite_master WHERE name="reaction"')

        if cur.fetchone()[0] == 0:  # no reaction table
            for init_command in init_commands:
                con.execute(init_command)  # Create tables
            con.commit()

        self.initialized = True

    def read(self, id, table='reaction'):
        """ Return an entire row of a table
        Parameters
        ---------

        id: int
            row integer
        table: str
            'reaction', 'publication', 'publication_system', 'reaction_system'
        """
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
        """
        Write publication info to db

        Parameters
        ----------
        values: dict with entries
            {'pub_id': str (short name for publication),
            'authors': list of str ()
            'journal': str,
            'volume': str,
            'number': str,
            'pages': 'str'
            'year': int,
            'publisher': str,
            'doi': str,
            'tags': list of str}
        """
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        values = (values['pub_id'],
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

        q = self.default + ',' + ', '.join('?' * len(values))
        cur.execute('INSERT OR IGNORE INTO publication VALUES ({})'.format(q),
                    values)

        pid = self.get_last_id(cur, table='publication')

        if self.connection is None:
            con.commit()
            con.close()

        return pid

    def write(self, values, data=None):
        """
        Write reaction info to db file

        Parameters
        ----------
        values: dict

        The values dict can include:
        {'chemical_composition': str (chemical composition on empty slab) ,
        'surface_composition': str (reduced chemical composition or
                                         shortname),
        'facet': str
        'sites': dict
            adsorption sites of species.
            f.ex: {'OH': 'ontop', 'O': 'hollow'}
        'coverages': dict
            coverage of adsorbates relative to the unit cell
            f.ex. {'OH': 0.25, 'O': 0.5})
        'reactants'/ 'products': dict
            keys with name of chemical species folloved by phase (gas, *)
            values are the prefactor in the reaction.
            For reaction H2Ogas -> 2Hstar + O star you would write:
            'reactants': {OHstar: 1, Hstar: 2}
            'products': {OHstar: 1, Hstar: 2}
        'reaction_energy': float
        'activation_energy': float
        'dft_code': str
        'dft_functional': str
        'username': str
        'pub_id': str
            Should match the pub_id of the corresponding publications
        }
        """
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        pub_id = values['pub_id']
        ase_ids = values['ase_ids']
        energy_corrections = values['energy_corrections']

        if ase_ids is not None:
            check_ase_ids(values, ase_ids)
        else:
            ase_ids = {}
        values = (values['chemical_composition'],
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

        """ Write to reaction table"""
        q = self.default + ',' + ', '.join('?' * len(values))
        cur.execute('INSERT INTO reaction VALUES ({})'.format(q),
                    values)
        id = self.get_last_id(cur)

        reaction_structure_values = []

        """ Write to publication_system and reaction_system tables"""
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

        return id

    def update(self, id, values, key_names='all'):
        """
        Update reaction info for a selected row

        Parameters
        ----------
        id: int
            row integer
        values: dict
            See write() method for details
        key_names: list or 'all'
            list with name of columns to update. Should match the keys-value
            pairs in values.
            default is 'all'
        """
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        pub_id = values['pub_id']
        ase_ids = values['ase_ids']
        energy_corrections = values['energy_corrections']

        if ase_ids is not None:
            check_ase_ids(values, ase_ids)
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

    def get_last_id(self, cur, table='reaction'):
        """
        Get the id of the last written row in table

        Parameters
        ----------
        cur: database connection().cursor() object
        table: str
            'reaction', 'publication', 'publication_system', 'reaction_system'

        Returns: id
        """
        cur.execute("SELECT seq FROM sqlite_sequence WHERE name='{0}'"
                    .format(table))
        result = cur.fetchone()
        if result is not None:
            id = result[0]
        else:
            id = 0
        return id

    def check(self, chemical_composition, reaction_energy):
        """
        Check if entry with same surface and energy is allready written
        to database file
        Parameters
        ----------
        chemcial_composition: str
        reaction_energy: str

        Returns id or None
        """
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
        """
        Check if entry with same surface and reaction is allready written
        to database file

        Parameters
        ----------
        chemcial_composition: str
        reactants: dict
        products: dict

        Returns id or None
        """
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

    def print_summary(self):
        self.stdout.write('------------------------------------------------\n')
        self.stdout.write('Reaction Summary: \n')
        self.stdout.write('------------------------------------------------\n')

        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        cur.execute("""
        SELECT
        surface_composition, reactants, products, reaction_energy,
        activation_energy, sites
        FROM
        reaction;""")
        rows = cur.fetchall()
        table = []
        for row in rows:
            equation = get_equation(json.loads(row[1]), json.loads(row[2]))
            table += [[row[0], equation, row[3], row[4], row[5]]]

        headers = ['Surface Composition', 'Equation', 'Reaction Energy',
                   'Activation Energy', 'Sites']
        self.stdout.write(tabulate(table, headers) + '\n')


def check_ase_ids(values, ase_ids):
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
                                      'ase_id', 'id'],
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
        if PY2:  # python 2
            if isinstance(v, unicode):
                v = v.encode('ascii', 'ignore')
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


def get_equation(reactants, products):
    equation = ''
    arrow = 0
    for column in (reactants, products):
        if arrow == 1:
            equation += ' -> '
        arrow += 1
        i = 0
        for key in sorted(column, key=len, reverse=True):
            prefactor = column[key]
            if 'gas' in key:
                key = key.replace('gas', '(g)')
            if 'star' in key:
                key = key.replace('star', '*')
            if not i == 0:
                if prefactor > 0:
                    equation += ' + '
                else:
                    equation += ' - '
                    prefactor *= -1
            if prefactor == 1:
                prefactor = ''

            equation += str(prefactor) + key
            i += 1
    return equation
