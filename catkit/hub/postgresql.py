import psycopg2
import os
import sys

try:
    from builtins import str as text
except BaseException:
    text = str

init_commands = [
    """CREATE TABLE publication (
    id SERIAL PRIMARY KEY,
    pub_id text UNIQUE,
    title text,
    authors jsonb,
    journal text,
    volume text,
    number text,
    pages text,
    year smallint,
    publisher text,
    doi text,
    tags jsonb
    );""",

    """CREATE TABLE publication_system (
    ase_id text REFERENCES systems(unique_id) ON DELETE CASCADE,
    pub_id text REFERENCES publication(pub_id) ON DELETE CASCADE,
    PRIMARY KEY (pub_id, ase_id)
    );""",

    """CREATE TABLE reaction (
    id SERIAL PRIMARY KEY,
    chemical_composition text,
    surface_composition text,
    facet text,
    sites jsonb,
    coverages jsonb,
    reactants jsonb,
    products jsonb,
    reaction_energy numeric,
    activation_energy numeric,
    dft_code text,
    dft_functional text,
    username text,
    pub_id text REFERENCES publication (pub_id) ON DELETE CASCADE
    );""",

    """CREATE TABLE reaction_system (
    name text,
    energy_correction numeric,
    ase_id text REFERENCES systems(unique_id) ON DELETE CASCADE,
    id integer REFERENCES reaction(id) ON DELETE CASCADE,
    PRIMARY KEY (id, ase_id)
    )"""
]

index_statements = [
    'CREATE INDEX idxpubid ON publication (pub_id);',
    'CREATE INDEX idxreacten ON reaction (reaction_energy);',
    'CREATE INDEX idxchemcomp ON reaction (chemical_composition);',
    'CREATE INDEX idxreact ON reaction USING GIN (reactants);',
    'CREATE INDEX idxprod ON reaction USING GIN (products);',
    'CREATE INDEX idxuser ON reaction (username);'
]

tsvector_statements = [
    """ALTER TABLE publication ADD COLUMN pubtextsearch tsvector;""",

    # Trigger doesn't work with all versions. Will make this work
    # later
    # """CREATE TRIGGER tsvectorupdatepub BEFORE INSERT OR UPDATE
    # ON publication FOR EACH ROW EXECUTE PROCEDURE
    # UPDATE publication SET pubtextsearch =
    # to_tsvector('english', coalesce(title, '') || ' ' ||
    # coalesce(authors::text, '') || ' ' || coalesce(year::text, '') || ' ' ||
    # coalesce(tags::text, ''))
    # ;""",
    # tsvector_update_trigger(pubtextsearch, 'pg_catalog.english',
    # title, authors, year, tags)

    """ALTER TABLE reaction ADD COLUMN textsearch tsvector;""",

    # """CREATE TRIGGER tsvectorupdate BEFORE INSERT OR UPDATE
    # ON reaction FOR EACH ROW EXECUTE PROCEDURE
    # tsvector_update_trigger(textsearch, 'pg_catalog.english',
    # chemical_compotision, facet, reactants, products);""",

    'CREATE INDEX idxsearch ON reaction USING GIN (textsearch);'
]
tsvector_update = [
    """UPDATE publication SET pubtextsearch =
    to_tsvector('simple', coalesce(title, '') || ' ' ||
    coalesce(authors::text, '') || ' ' || coalesce(year::text, '') || ' ' ||
    coalesce(tags::text, ''));
    """,

    """
    UPDATE reaction SET textsearch =
    to_tsvector('simple', coalesce(regexp_replace(
    regexp_replace(chemical_composition,
    '([0-9])', '', 'g'), '()([A-Z])', '\1 \2','g'), '') || ' ' ||
    coalesce(facet, '') || ' ' ||
    replace(replace(coalesce(reactants::text, '') || ' ' ||
    coalesce(products::text, ''), 'star',''), 'gas', ''));
    """
]


class CathubPostgreSQL:
    """ Class for setting up the catalysis hub reaction energy database
    on postgreSQL server.
    """

    def __init__(self, user='catroot', password=None, stdin=sys.stdin,
                 stdout=sys.stdout):
        self.initialized = False
        self.connection = None
        self.id = None
        if user == 'catroot' or user == 'catvisitor':
            self.schema = 'public'
        else:
            self.schema = user
        self.user = user
        self.server = 'catalysishub.c8gwuc8jwb7l.us-west-2.rds.amazonaws.com'
        if password is None:
            password = os.environ['DB_PASSWORD']
        self.password = password
        self.stdin = stdin
        self.stdout = stdout

    def _connect(self):
        con = psycopg2.connect(host=self.server,
                               user=self.user,
                               password=self.password,
                               port=5432,
                               database='catalysishub')

        return con

    def __enter__(self):
        assert self.connection is None
        self.connection = self._connect()
        return self

    def __exit__(self, exc_type):  # , exc_value, tb):
        if exc_type is None:
            self.connection.commit()
        else:
            self.connection.rollback()
        self.connection.close()
        self.connection = None

    def _initialize(self, con):
        if self.initialized:
            return
        cur = con.cursor()

        self.stdout.write("_initialize start\n")

        set_schema = 'SET search_path = {0};'.format(self.schema)
        cur.execute(set_schema)

        self.stdout.write(
            "_initialize set schema to {self.schema}\n".format(**locals()))

        from ase.db.postgresql import PostgreSQLDatabase
        PostgreSQLDatabase()._initialize(con)

        cur.execute("""SELECT to_regclass('publication');""")
        if cur.fetchone()[0] is None:  # publication doesn't exist
            for init_command in init_commands:
                self.stdout.write(init_command + '\n')
                cur.execute(init_command)
            for statement in index_statements:
                self.stdout.write(statement + '\n')
                cur.execute(statement)
            for statement in tsvector_statements:
                self.stdout.write(statement + '\n')
                cur.execute(statement)
            con.commit()
        self.initialized = True
        return self

    def create_user(self, user):
        from pwgen import pwgen
        con = self.connection or self._connect()
        cur = con.cursor()

        cur.execute('CREATE SCHEMA {0};'.format(user))
        # self._initialize(schema=schema_name)
        password = pwgen(8)
        cur.execute(
            "CREATE USER {0} with PASSWORD '{1}';".format(user, password))
        cur.execute(
            'GRANT ALL PRIVILEGES ON SCHEMA {0} TO {1};'.format(user))
        cur.execute('GRANT USAGE ON SCHEMA public TO {0};'.format(user))
        cur.execute(
            'GRANT SELECT ON ALL TABLES IN SCHEMA public TO {0};'.format(user))
        cur.execute(
            'ALTER ROLE {0} SET search_path TO {0};'.format(user))

        self.stdout.write(
            'CREATED USER {0} WITH PASSWORD {1}\n'.format(user, password))

        con.commit()
        con.close()

        self.schema = user
        self.user = user
        self.password = password
        con = self._connect()
        self._initialize(con)

        con.commit()
        con.close()

        return self

    def remove_user(self, user):
        assert self.user == 'catroot'
        con = self.connection or self._connect()
        cur = con.cursor()
        cur.execute('DROP SCHEMA {0} CASCADE;'.format(user))
        cur.execute('REVOKE USAGE ON SCHEMA public FROM {0};'.format(user))
        cur.execute(
            'REVOKE SELECT ON ALL TABLES IN SCHEMA public FROM {0};'
            .format(user))
        cur.execute(
            'DROP ROLE {0};'.format(user))
        self.stdout.write(
            'REMOVED USER {0}\n'.format(user))
        con.commit()
        con.close()

        return self

    def status(self, table='reaction'):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        cur.execute("SELECT COUNT(id) from {0};".format(table))
        count = cur.fetchone()
        return count[0]

    def read(self, id, table='reaction'):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        cur.execute(
            """SELECT column_name FROM information_schema.columns
            WHERE table_schema = 'public' AND table_name='{0}';"""
            .format(table))
        columns = cur.fetchall()

        if id == 'all':
            cur.execute('SELECT * FROM \n {0} \n'.format(table,
                                                         table))
        else:
            cur.execute('SELECT * FROM \n {0} \n WHERE \n {1}.id={2}'
                        .format(table, table, id))
        row = cur.fetchall()

        return columns, row

    def write_publication(self, pub_values):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        # pub_id = pub_values[1].encode('ascii','ignore')
        pub_id = pub_values[1]
        cur.execute(
            """SELECT id from publication where pub_id='{0}'""".format(pub_id))
        row = cur.fetchone()
        if row is not None:  # len(row) > 0:
            id = row  # [0]
        else:
            key_str, value_str = get_key_value_str(pub_values, 'publication')
            insert_command = """INSERT INTO publication ({0}) VALUES
            ({1}) RETURNING id;""".format(key_str, value_str)

            cur.execute(insert_command)
            id = cur.fetchone()[0]

        if self.connection is None:
            con.commit()
            con.close()
        return id, pub_id

    def write(self, values, table='reaction'):

        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        key_str, value_str = get_key_value_str(values, table)

        insert_command = 'INSERT INTO {0} ({1}) VALUES ({2}) RETURNING id;'\
            .format(table, key_str, value_str)

        cur.execute(insert_command)
        id = cur.fetchone()[0]

        if self.connection is None:
            con.commit()
            con.close()
        return id

    def update(self, id, values, key_names='all'):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        key_str, value_str = get_key_value_str(values)

        update_command = 'UPDATE reaction SET ({0}) = ({1}) WHERE id = {2};'\
            .format(key_str, value_str, id)

        cur.execute(update_command)

        if self.connection is None:
            con.commit()
            con.close()
        return id

    def update_publication(self, pub_dict):
        import json
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        pub_id = pub_dict['pub_id']

        values = pub_dict.values()
        key_str = ', '.join(pub_dict.keys())
        value_str = "'{0}'".format(values[0])
        for v in values[1:]:
            if isinstance(v, text):
                v = v.encode('utf8', 'ignore')
            if v is None or v == '':
                value_str += ", {0}".format('NULL')
            elif isinstance(v, str):
                value_str += ", '{0}'".format(v)
            elif isinstance(v, list):
                value_str += ", '{0}'".format(json.dumps(v))
            else:
                value_str += ", {0}".format(v)

        update_command = \
            """UPDATE publication SET ({0}) = ({1}) WHERE pub_id='{2}';"""\
            .format(key_str, value_str, pub_id)

        self.stdout.write(update_command + '\n')
        cur.execute(update_command)

        if self.connection is None:
            con.commit()
            con.close()

        return

    def delete(self, authorlist, year, doi=None):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        if doi is None:
            delete_command =\
                """DELETE from reaction
                WHERE
                publication -> 'authors' = '{0}' and year = {2};"""\
                    .format(authorlist, year)
        cur.execute(delete_command)
        count = cur.fetchone()[0]
        if self.connection is None:
            con.commit()
            con.close()

        return count

    def transfer(self, filename_sqlite, start_id=1, write_ase=True,
                 write_publication=True, write_reaction=True,
                 write_reaction_system=True, block_size=1000,
                 start_block=0):

        self.stdout.write('Starting transfer\n')
        con = self.connection or self._connect()
        self._initialize(con)
        self.stdout.write('Finished initialization\n')
        cur = con.cursor()
        self.stdout.write('Got a cursor\n')

        set_schema = 'SET search_path = {0};'.format(self.schema)
        cur.execute(set_schema)

        import os
        import time
        self.stdout.write('Imported os\n')

        import ase.db
        self.stdout.write('Imported ase.db\n')
        self.stdout.write('Building server_name\n')
        server_name = "postgres://{0}:{1}@{2}:5432/catalysishub".format(
            self.user, self.password, self.server)

        self.stdout.write('Connecting to {server_name}\n'.format(**locals()))

        nrows = 0
        if write_ase:
            print('Transfering atomic structures')
            db = ase.db.connect(filename_sqlite)
            n_structures = db.count()
            n_blocks = int(n_structures / block_size) + 1
            t_av = 0
            for block_id in range(start_block, n_blocks):
                i = block_id - start_block
                t1 = time.time()
                b0 = block_id * block_size + 1
                b1 = (block_id + 1) * block_size + 1
                # self.stdout.write(str(block_id) + ' ' +
                #                  'from ' + str(b0) + ' to ' + str(b1) + '\n')
                if block_id + 1 == n_blocks:
                    b1 = n_structures + 1

                rows = list(db.select('{}<id<{}'.format(b0 - 1, b1)))

                with ase.db.connect(server_name, type='postgresql') as db2:

                    db2.write(rows)

                nrows += len(rows)
                t2 = time.time()
                dt = t2 - t1
                t_av = (t_av * i + dt) / (i + 1)

                self.stdout.write(
                    '  Finnished Block {0} / {1} in {2} sec'
                    .format(block_id + 1, n_blocks, dt))
                self.stdout.write(
                    '    Completed transfer of {0} atomic structures.'
                    .format(nrows))
                self.stdout.write('    Estimated time left: {0} sec'.format(
                    t_av * (n_blocks - block_id)))

        from catkit.hub.cathubsqlite import CathubSQLite
        db = CathubSQLite(filename_sqlite)
        con_lite = db._connect()
        cur_lite = con_lite.cursor()

        # write publication
        Npub = 0
        Npubstruc = 0
        if write_publication:
            try:
                npub = db.get_last_pub_id(cur_lite)
            except BaseException:
                npub = 1
            for id_lite in range(1, npub + 1):
                Npub += 1
                row = db.read(id=id_lite, table='publication')
                if len(row) == 0:
                    continue
                values = row[0]
                pid, pub_id = self.write_publication(values)

            # Publication structures connection
            cur_lite.execute("""SELECT * from publication_system;""")
            rows = cur_lite.fetchall()
            for row in rows:
                Npubstruc += 1
                values = row[:]
                key_str, value_str = get_key_value_str(
                    values, table='publication_system')

                set_schema = 'SET search_path = {0};'.format(self.schema)
                cur.execute(set_schema)
                print("[SET SCHEMA] {set_schema}".format(**locals()))

                insert_command = """INSERT INTO publication_system ({0})
                VALUES ({1}) ON CONFLICT DO NOTHING;"""\
                    .format(key_str, value_str)

                cur.execute(insert_command)
                # self.write(values, table='publication_system')
            con.commit()

        Ncat = 0
        Ncatstruc = 0

        if write_reaction:
            n = db.get_last_id(cur_lite)
            select_ase = """SELECT * from reaction_system where id={};"""
            for id_lite in range(start_id, n + 1):
                row = db.read(id_lite)
                if len(row) == 0:
                    continue
                values = row[0]

                id = self.check(values[13], values[1], values[6], values[7],
                                values[8], strict=True)
                update_rs = False

                if id is not None:
                    id = self.update(id, values)
                    self.stdout.write(
                        'Updated reaction db with row id = {}\n'.format(id))
                    update_rs = True
                else:
                    Ncat += 1
                    id = self.write(values)
                    self.stdout.write(
                        'Written to reaction db row id = {0}\n'.format(id))

                cur_lite.execute(select_ase.format(id_lite))
                rows = cur_lite.fetchall()
                if write_reaction_system:
                    if update_rs:
                        cur.execute("""Delete from reaction_system231
                        where reaction_id={0}""".format(id))
                    for row in rows:
                        Ncatstruc += 1
                        values = list(row)
                        if len(values) == 3:
                            values.insert(1, None)

                        values[3] = id

                        key_str, value_str = \
                            get_key_value_str(values, table='reaction_system')

                        set_schema = 'SET search_path = {0};'.format(
                            self.schema)
                        cur.execute(set_schema)
                        # print("[SET SCHEMA] {set_schema}".format(**locals()))

                        insert_command = """INSERT INTO reaction_system
                        ({0}) VALUES ({1}) ON CONFLICT DO NOTHING;"""\
                            .format(key_str, value_str)

                        # print(
                        #     "[INSERT COMMAND] {insert_command}"
                        #     .format(**locals()))
                        cur.execute(insert_command)

                con.commit()  # Commit reaction_system for each row

        for statement in tsvector_update:
            cur.execute(statement)

        if self.connection is None:
            con.commit()
            con.close()

        self.stdout.write('Inserted into:\n')
        self.stdout.write('  systems: {0}\n'.format(nrows))
        self.stdout.write('  publication: {0}\n'.format(Npub))
        self.stdout.write('  publication_system: {0}\n'.format(Npubstruc))
        self.stdout.write('  reaction: {0}\n'.format(Ncat))
        self.stdout.write('  reaction_system: {0}\n'.format(Ncatstruc))

    def check(self, pub_id, chemical_composition, reactants, products,
              reaction_energy=None, strict=True):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        keys = 'pub_id, chemical_composition,  reactants, products'
        values = [pub_id, chemical_composition, reactants, products]
        placeholder = """'{}', '{}', '{}', '{}'"""
        if strict:
            assert reaction_energy is not None
            placeholder += ", {}"
            keys += ', reaction_energy'
            values.append(reaction_energy)

        placeholder += """);"""
        arguments = [keys] + values

        statement = \
            """SELECT id
        FROM reaction WHERE
        ({}) =
        (""" + placeholder

        statement = statement.format(*arguments)

        cur.execute(statement)
        rows = cur.fetchall()

        if len(rows) > 0:
            id = rows[0][0]
        else:
            id = None
        return id

    def publication_status(self):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        select_statement = \
            """SELECT distinct publication
        FROM reaction WHERE
        publication ->> 'doi' is null
        OR
        publication -> 'doi' is null;"""
        cur.execute(select_statement)
        pubs = cur.fetchall()

        return pubs


def get_key_value_str(values, table='reaction'):
    key_str = {'reaction': """chemical_composition, surface_composition,
    facet, sites, coverages, reactants, products, reaction_energy,
    activation_energy, dft_code, dft_functional, username, pub_id""",
               'publication': """pub_id, title, authors, journal,
               volume, number, pages, year, publisher, doi, tags""",
               'reaction_system': """name, energy_correction, ase_id, id""",
               'publication_system': """ase_id, pub_id"""}

    start_index = 1
    if table == 'publication_system' or table == 'reaction_system':
        start_index = 0
    value_str = "'{0}'".format(values[start_index])
    for v in values[start_index + 1:]:
        # print("\n\n\nDIR TYPE {v}".format(**locals()))
        # print(dir(v))
        # print(type(v))
        # if isinstance(v, text):
        # v = v.encode('utf8','ignore')
        # print("ISINSTANCE TEXT {v}".format(**locals()))
        # elif hasattr(v, 'encode'):
        # v = v.encode('utf8','ignore')
        # print("HASATTR ENCODE {v}".format(**locals()))

        if v is None or v == '':
            value_str += ", {0}".format('NULL')
        elif isinstance(v, str):
            value_str += ", '{0}'".format(v)
        else:
            value_str += ", {0}".format(v)
        # print(value_str)

    return key_str[table], value_str
