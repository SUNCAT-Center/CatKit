import os
import sys
import time
import json
import random
import psycopg2
from psycopg2.extras import execute_values
import ase.db
from ase.db.core import now
from ase.db.postgresql import PostgreSQLDatabase
from past.utils import PY2

from .cathubsqlite import CathubSQLite

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


def pwgen(n):
    return ''.join([
        random.choice(
            'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            'abcdefghijklmnopqrstuvwxyz'
            '0123456789'
        ) for i in range(n)])


class CathubPostgreSQL:
    """ Class for setting up the catalysis hub reaction energy database
    on postgreSQL server.
    """

    def __init__(self, user='catroot', password=None, stdin=sys.stdin,
                 stdout=sys.stdout):
        self.initialized = False
        self.connection = None
        self.id = None
        self.server = 'catalysishub.c8gwuc8jwb7l.us-west-2.rds.amazonaws.com'
        self.database = 'catalysishub'
        if user == 'catroot' or user == 'catvisitor':
            self.schema = 'public'
            if password is None:
                password = os.environ['DB_PASSWORD']
        elif user == 'postgres':  # For testing on travis
            self.schema = 'public'
            self.server = 'localhost'
            self.database = 'travis_ci_test'
            self.password = ''
        else:
            self.schema = user
        self.user = user
        self.password = password
        self.stdin = stdin
        self.stdout = stdout

    def _connect(self):
        con = psycopg2.connect(host=self.server,
                               user=self.user,
                               password=self.password,
                               port=5432,
                               database=self.database)
        self.server_name = "postgres://{0}:{1}@{2}:5432/{3}".format(
            self.user, self.password, self.server, self.database)

        return con

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
        cur = con.cursor()

        self.stdout.write("_initialize start\n")

        set_schema = 'SET search_path TO {0};'\
                     .format(self.schema)
        cur.execute(set_schema)
        con.commit()

        cur.execute("show search_path;")
        schema = cur.fetchone()[0].split(', ')

        self.stdout.write(
            "_initialize set schema to {self.schema}\n".format(**locals()))

        PostgreSQLDatabase()._initialize(con)

        cur.execute("""SELECT to_regclass('publication');""")
        if cur.fetchone()[0] is None:  # publication doesn't exist
            self.stdout.write("_initialize create tables:\n")
            for init_command in init_commands:
                self.stdout.write(init_command + '\n')
                cur.execute(init_command)
            self.stdout.write("_initialize create indexes:\n")
            for statement in index_statements:
                self.stdout.write(statement + '\n')
                cur.execute(statement)
            self.stdout.write("_initialize create text search columns:\n")
            for statement in tsvector_statements:
                self.stdout.write(statement + '\n')
                cur.execute(statement)
            con.commit()
        self.initialized = True
        return self

    def get_ase_db(self):
        if not self.connection:
            self._connect()
        print(self.server_name)
        return ase.db.connect(self.server_name)

    def create_user(self, user, table_privileges=['ALL PRIVILEGES'],
                    schema_privileges=['ALL PRIVILEGES'],
                    row_limit=50000):
        con = self.connection or self._connect()
        cur = con.cursor()

        cur.execute('CREATE SCHEMA {0};'.format(user))
        # self._initialize(schema=schema_name)
        password = pwgen(8)
        cur.execute(
            "CREATE USER {user} with PASSWORD '{password}';"
            .format(user=user, password=password))

        """ Grant SELECT on public schema """
        cur.execute('GRANT USAGE ON SCHEMA public TO {user};'
                    .format(user=user))
        cur.execute(
            'GRANT SELECT ON ALL TABLES IN SCHEMA public TO {user};'
            .format(user=user))
        cur.execute(
            'ALTER ROLE {user} SET search_path TO {user};'
            .format(user=user))

        self.stdout.write(
            'CREATED USER {user} WITH PASSWORD {password}\n'
            .format(user=user, password=password))

        """ initialize user-schema """
        old_schema = self.schema
        self.initialized = False
        self.schema = user
        self._initialize(con)

        """ Privileges on user-schema"""
        cur.execute(
            'GRANT {privileges} ON SCHEMA {user} TO {user};'
            .format(privileges=', '.join(schema_privileges), user=user))
        cur.execute(
            'GRANT {privileges} ON ALL TABLES IN SCHEMA {user} TO {user};'
            .format(privileges=', '.join(table_privileges), user=user))
        cur.execute(
            'GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA {user} TO {user};'
            .format(user=user))
        con.commit()

        if row_limit:
            """ Limit number of rows"""
            for table in ['reaction', 'publication', 'systems',
                          'reaction_system', 'publication_system',
                          'information']:
                table_factor = 1
                if table in ['reaction_system', 'publication_system']:
                    table_factor = 15
                elif table == 'publication':
                    table_factor = 1 / 100
                elif table == 'information':
                    table_factor = 1 / 100

                trigger_function = """
                CREATE OR REPLACE FUNCTION
                check_number_of_rows_{user}_{table}()
                RETURNS TRIGGER AS
                $BODY$
                BEGIN
                    IF (SELECT count(*) FROM {user}.{table}) > {row_limit}
                    THEN
                        RAISE EXCEPTION
                            'INSERT statement exceeding maximum number of rows';
                    END IF;
                    RETURN NEW;
                END;
                $BODY$
                LANGUAGE plpgsql""".format(user=user, table=table,
                                           row_limit=row_limit * table_factor)
                cur.execute(trigger_function)

                trigger = """
                DROP TRIGGER IF EXISTS tr_check_number_of_rows_{user}_{table}
                    on {user}.{table};
                CREATE TRIGGER tr_check_number_of_rows_{user}_{table}
                BEFORE INSERT ON {user}.systems
                FOR EACH ROW EXECUTE PROCEDURE check_number_of_rows_{user}_{table}();
                """.format(user=user, table=table)
                cur.execute(trigger)

        self.schema = old_schema
        set_schema = 'ALTER ROLE {user} SET search_path TO {schema};'\
                     .format(user=self.user, schema=self.schema)
        cur.execute(set_schema)

        if self.connection is None:
            con.commit()
            con.close()

        return password

    def delete_user(self, user):
        """ Delete user and all data"""
        assert self.user == 'catroot' or self.user == 'postgres'
        assert not user == 'public'
        con = self.connection or self._connect()
        cur = con.cursor()
        cur.execute('DROP SCHEMA {user} CASCADE;'.format(user=user))
        cur.execute('REVOKE USAGE ON SCHEMA public FROM {user};'
                    .format(user=user))
        cur.execute(
            'REVOKE SELECT ON ALL TABLES IN SCHEMA public FROM {user};'
            .format(user=user))
        cur.execute(
            'DROP ROLE {user};'.format(user=user))
        self.stdout.write(
            'REMOVED USER {user}\n'.format(user=user))

        if self.connection is None:
            con.commit()
            con.close()

        return self

    def release(self, pub_ids=None, userhandle=None, from_schema='upload',
                to_schema='public'):
        """ Transfer dataset from one schema to another"""

        assert pub_ids or userhandle,\
            "Specify either pub_ids or userhandle"
        assert not (pub_ids and userhandle),\
            "Specify either pub_ids or userhandle"

        con = self.connection or self._connect()
        cur = con.cursor()
        assert self.user in ['release', 'catroot', 'postgres'], \
            "You don't have permission to perform this operation"

        if userhandle:
            cur.execute(
                """SELECT distinct pub_id
                FROM {from_schema}.reaction
                WHERE username = '{username}'"""
                .format(from_schema=from_schema,
                        username=userhandle))

            pub_ids = [id[0] for id in cur.fetchall()]

        for pub_id in pub_ids:
            self.stdout.write(
                """Releasing publication {pub_id} from
                {from_schema} to {schema} \n"""
                .format(pub_id=pub_id, from_schema=from_schema,
                        schema=to_schema))

            mtime = now()
            cur.execute(
                """UPDATE {from_schema}.systems SET
                mtime = {mtime}
                WHERE unique_id in
                  (SELECT distinct ase_id
                   FROM {from_schema}.publication_system
                   WHERE pub_id = '{pub_id}')"""
                .format(from_schema=from_schema,
                        mtime=mtime, pub_id=pub_id))

            columns = get_key_str('systems', start_index=1)
            cur.execute(
                """INSERT INTO {schema}.systems ({columns})
                SELECT {columns}
                FROM {from_schema}.systems
                WHERE unique_id in
                  (SELECT distinct ase_id
                   FROM {from_schema}.publication_system
                   WHERE pub_id = '{pub_id}')"""
                .format(from_schema=from_schema,
                        schema=to_schema,
                        columns=columns,
                        pub_id=pub_id))

            columns = get_key_str('publication', start_index=1)  # new id
            cur.execute(
                """INSERT INTO {schema}.publication ({columns})
                SELECT {columns}
                FROM {from_schema}.publication
                WHERE pub_id = '{pub_id}'"""
                .format(from_schema=from_schema,
                        schema=to_schema, columns=columns,
                        pub_id=pub_id))
            cur.execute(
                """INSERT INTO {schema}.publication_system
                SELECT *
                FROM {from_schema}.publication_system
                WHERE pub_id = '{pub_id}'"""
                .format(from_schema=from_schema,
                        schema=to_schema, columns=columns, pub_id=pub_id))

            columns = get_key_str('reaction', start_index=1)  # new id

            cur.execute(
                """INSERT INTO {schema}.reaction ({columns})
                SELECT {columns}
                FROM {from_schema}.reaction
                WHERE pub_id = '{pub_id}'
                ORDER BY {from_schema}.reaction.id
                RETURNING id"""
                .format(from_schema=from_schema,
                        schema=to_schema, columns=columns, pub_id=pub_id))

            new_ids = [id[0] for id in cur.fetchall()]

            cur.execute(
                """SELECT * from {from_schema}.reaction_system
                WHERE ase_id in
                (SELECT distinct ase_id
                FROM {from_schema}.publication_system
                WHERE pub_id = '{pub_id}') ORDER BY id"""
                .format(from_schema=from_schema,
                        pub_id=pub_id))

            reaction_system0 = cur.fetchall()
            reaction_system_values = []
            id0 = reaction_system0[0][3]
            i = 0
            for row in reaction_system0:
                row = list(row)
                if not id0 == row[3]:
                    i += 1
                    id0 = row[3]
                row[3] = new_ids[i]
                reaction_system_values += [tuple(row)]
            key_str = get_key_str('reaction_system')
            insert_command = """
            INSERT INTO {schema}.reaction_system ({key_str})
            VALUES %s ON CONFLICT DO NOTHING;"""\
            .format(schema=to_schema, key_str=key_str)

            execute_values(cur=cur, sql=insert_command,
                           argslist=reaction_system_values, page_size=1000)
            self.stdout.write('Transfer complete\n')

        if self.user == 'catroot':
            self.delete_publication(pub_id, schema='upload')

        if self.connection is None:
            con.commit()
            con.close()

        return

    def get_pub_id_owner(self, pub_id):
        """Check if a user owns a publication"""
        con = self.connection or self._connect()
        cur = con.cursor()

        cur.execute(
            """SELECT distinct username from upload.reaction
            WHERE pub_id = '{pub_id}'"""
            .format(pub_id=pub_id))

        username = cur.fetchall()[0][0]

        return username

    def delete_publication(self, pub_id, schema='upload'):
        """ Delete dataset from upload schema"""
        if schema == 'upload':
            user = 'upload_admin'
        elif schema == 'public':
            user = 'catroot'

        if not self.user == 'catroot':
            assert self.user == user, \
                "You don't have permission to perform this operation"

        con = self.connection or self._connect()
        cur = con.cursor()

        self.stdout.write('Deleting publication: {pub_id} from {schema}\n'
                          .format(pub_id=pub_id, schema=schema))

        cur.execute("""SELECT to_regclass('keys');""")
        if cur.fetchone()[0] is not None:  # remove data from old tables
            old_tables = ['text_key_values', 'number_key_values',
                          'species', 'keys']
            for table in old_tables:
                cur.execute(
                    """DELETE FROM {schema}.{table}"""
                    .format(schema=schema,
                            table=table))

        cur.execute(
            """DELETE FROM {schema}.systems
            WHERE unique_id in
            (SELECT distinct ase_id
            FROM {schema}.publication_system
            WHERE pub_id = '{pub_id}')"""
            .format(schema=schema,
                    pub_id=pub_id))

        cur.execute(
            """ DELETE FROM {schema}.publication
            WHERE pub_id = '{pub_id}'"""
            .format(schema=schema,
                    pub_id=pub_id))

        self.stdout.write('Delete complete\n')

        if self.connection is None:
            con.commit()
            con.close()
        return

    def status(self, table='reaction'):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        cur.execute("SELECT COUNT(id) from {table};".format(table=table))
        count = cur.fetchone()
        return count[0]

    def read(self, id, table='reaction'):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        cur.execute(
            """SELECT column_name FROM information_schema.columns
            WHERE table_schema = 'public' AND table_name='{table}';"""
            .format(table=table))
        columns = cur.fetchall()

        if id == 'all':
            cur.execute('SELECT * FROM \n {table} \n'.format(table=table))
        else:
            cur.execute('SELECT * FROM \n {table} \n WHERE \n {table}.id={id}'
                        .format(table=table, id=id))
        row = cur.fetchall()

        return columns, row

    def write_publication(self, pub_values):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        # pub_id = pub_values[1].encode('ascii','ignore')
        pub_id = pub_values[1]
        cur.execute(
            """SELECT id from publication where pub_id='{pub_id}'"""
            .format(pub_id=pub_id))
        row = cur.fetchone()
        if row is not None:  # len(row) > 0:
            id = row  # [0]
        else:
            key_str = get_key_str('publication', start_index=1)
            value_str = get_value_str(pub_values, start_index=1)

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

        key_str = get_key_str(table, start_index=1)
        value_str = get_value_str(values, start_index=1)

        insert_command = 'INSERT INTO {0} ({1}) VALUES ({2}) RETURNING id;'\
            .format(table, key_str, value_str)

        cur.execute(insert_command)
        id = cur.fetchone()[0]

        if self.connection is None:
            con.commit()
            con.close()
        return id

    def write_reaction(self, value_dict):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        ase_ids = value_dict['ase_ids']
        energy_corrections = value_dict.get('energy_corrections', {})

        key_list = get_key_list(start_index=1)
        values = [value_dict[key] for key in key_list]

        key_str = get_key_str('reaction', start_index=1)
        value_str = get_value_str(values)

        insert_command = \
            """INSERT INTO reaction ({0})
            VALUES ({1}) RETURNING id;"""\
            .format(key_str, value_str)

        cur.execute(insert_command)
        id = cur.fetchone()[0]

        reaction_system_values = []

        """ Write to reaction_system tables"""
        for name, ase_id in ase_ids.items():
            if name in energy_corrections:
                energy_correction = energy_corrections[name]
            else:
                energy_correction = 0

            reaction_system_values += [tuple([name, energy_correction,
                                              ase_id, id])]

        key_str = get_key_str('reaction_system')
        insert_command = """INSERT INTO reaction_system
        ({0}) VALUES %s ON CONFLICT DO NOTHING;""".format(key_str)

        execute_values(cur=cur, sql=insert_command,
                       argslist=reaction_system_values, page_size=1000)

        if self.connection is None:
            con.commit()
            con.close()

        return id

    def update_reaction(self, id, ase_ids=None, energy_corrections={},
                        **kwargs):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        key_str = ', '.join(list(kwargs.keys()))
        value_str = get_value_str(list(kwargs.values()), start_index=0)

        update_command = 'UPDATE reaction SET ({0}) = ({1}) WHERE id = {2};'\
            .format(key_str, value_str, id)

        cur.execute(update_command)

        if ase_ids:
            delete_command = 'DELETE from reaction_system where id = {};'\
                .format(id)
            cur.execute(delete_command)

            """ Write to reaction_system tables"""
            reaction_system_values = []
            for name, ase_id in ase_ids.items():
                if name in energy_corrections:
                    energy_correction = energy_corrections[name]
                else:
                    energy_correction = 0

                    reaction_system_values += [tuple([name, energy_correction,
                                                      ase_id, id])]

            key_str = get_key_str('reaction_system')
            insert_command = """INSERT INTO reaction_system
            ({0}) VALUES %s ON CONFLICT DO NOTHING;""".format(key_str)

            execute_values(cur=cur, sql=insert_command,
                           argslist=reaction_system_values, page_size=1000)

        if self.connection is None:
            con.commit()
            con.close()
        return id

    def delete_reaction(self, id):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        delete_command = 'DELETE FROM reaction where id = {};'.format(id)
        cur.execute(delete_command)

        delete_command = 'DELETE from reaction_system where id = {};'\
            .format(id)
        cur.execute(delete_command)

        if self.connection is None:
            con.commit()
            con.close()
        return id

    def update_publication(self, pub_dict):
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

    def truncate_schema(self):
        """ Will delete all data in schema. Only for test use!"""

        assert self.server == 'localhost'
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()

        cur.execute('DELETE FROM publication;')
        cur.execute('TRUNCATE systems CASCADE;')

        con.commit()
        con.close()

        return

    def transfer(self, filename_sqlite, block_size=1000,
                 start_block=0, write_ase=True,
                 write_publication=True, write_reaction=True,
                 write_reaction_system=True, check=False):
        """ Transfer data from local sqlite3 .db file to the
        catalysis-hub postgreSQL server

        Parameters:
        filename_sqlite: str
            name of .db file
        block_size: int (default 1000)
            Number of atomic structures and reactions to write together
            in each block.
        start_block: int (default 0)
            Block to start with
        write_ase: bool
            whether or not to write atomic structures
        write_publication: bool
            whether or not to transfer publication table
        write_reaction: bool
            whether or not to transfer reaction table
        write_reaction_system: bool
            whether or not to write reaction_system table
        """

        self.stdout.write('Starting transfer\n')
        con = self.connection or self._connect()
        self._initialize(con)
        self.stdout.write('Finished initialization\n')
        cur = con.cursor()
        self.stdout.write('Got a cursor\n')
        self.stdout.write('Connecting to {0}\n'.format(self.server_name))

        nrows = 0
        if write_ase:
            self.stdout.write('Transfering atomic structures\n')
            db = ase.db.connect(filename_sqlite)
            n_structures = db.count()
            n_blocks = n_structures // block_size + 1
            t_av = 0
            for block_id in range(start_block, n_blocks):
                i = block_id - start_block
                t1 = time.time()
                b0 = block_id * block_size
                b1 = (block_id + 1) * block_size + 1

                if block_id + 1 == n_blocks:
                    b1 = n_structures + 1

                rows = list(db.select('{}<id<{}'.format(b0, b1)))

                with ase.db.connect(self.server_name, type='postgresql') as db2:
                    # write one row at the time until ase is updated
                    # db2.write(rows)
                    for row in rows:
                        db2.write(row)

                nrows += len(rows)
                t2 = time.time()
                dt = t2 - t1
                t_av = (t_av * i + dt) / (i + 1)

                self.stdout.write(
                    '  Finnished Block {0} / {1} in {2} sec\n'
                    .format(block_id + 1, n_blocks, dt))
                self.stdout.write(
                    '    Completed transfer of {0} atomic structures\n'
                    .format(nrows))
                self.stdout.write('    Estimated time left: {0} sec\n'.format(
                    t_av * (n_blocks - block_id - 1)))

        db = CathubSQLite(filename_sqlite)
        con_lite = db._connect()
        cur_lite = con_lite.cursor()

        Npub = 0
        Npubstruc = 0
        if write_publication:
            self.stdout.write('Transfering publications\n')
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
            publication_system_values = []
            rows = cur_lite.fetchall()
            for row in rows:
                Npubstruc += 1
                values = list(row)
                value_list = get_value_list(values)
                publication_system_values += [tuple(value_list)]

            # Insert into publication_system table
            key_str = get_key_str(table='publication_system')
            insert_command = """INSERT INTO publication_system ({0})
            VALUES %s ON CONFLICT DO NOTHING;"""\
                .format(key_str)

            execute_values(cur=cur, sql=insert_command,
                           argslist=publication_system_values, page_size=1000)

            # Write pub_id to systems table
            cur.execute("""UPDATE systems SET
            key_value_pairs=jsonb_set(key_value_pairs, '{{"pub_id"}}', '"{pub_id}"')
            WHERE unique_id IN
            (SELECT ase_id from publication_system WHERE pub_id='{pub_id}')"""
                        .format(pub_id=pub_id))

            con.commit()
            self.stdout.write('  Completed transfer of publications\n')

        Ncat = 0
        Ncatstruc = 0

        if write_reaction:
            self.stdout.write('Transfering reactions')
            cur.execute('SELECT max(id) from reaction;')
            ID = cur.fetchone()[0] or 0

            n_react = db.get_last_id(cur_lite)

            n_blocks = int(n_react / block_size) + 1
            t_av = 0
            for block_id in range(start_block, n_blocks):
                reaction_values = []
                reaction_system_values = []
                Ncat0 = Ncat
                Ncatstruc0 = Ncatstruc

                i = block_id - start_block
                t1 = time.time()
                b0 = block_id * block_size + 1
                b1 = (block_id + 1) * block_size + 1
                if block_id + 1 == n_blocks:
                    b1 = n_react + 1

                for id_lite in range(b0, b1):
                    row = db.read(id_lite)
                    if len(row) == 0:
                        continue
                    values = row[0]

                    id = None
                    update_rs = False
                    if id is not None:
                        id = self.update(id, values)
                        self.stdout.write(
                            'Updated reaction db with row id = {}\n'
                            .format(id))
                        update_rs = True
                    else:
                        ID += 1
                        Ncat += 1
                        value_list = get_value_list(values)
                        value_list[0] = ID  # set new ID
                        reaction_values += [tuple(value_list)]
                        if write_reaction_system:
                            cur_lite.execute("SELECT * from reaction_system where id={};"
                                             .format(id_lite))
                            rows = cur_lite.fetchall()
                            if update_rs:
                                cur.execute("""Delete from reaction_system
                                where id={0}""".format(id))
                            for row in rows:
                                Ncatstruc += 1
                                values = list(row)
                                if len(values) == 3:
                                    values.insert(1, None)
                                value_list = get_value_list(values)
                                value_list[3] = ID
                                reaction_system_values += [tuple(value_list)]

                q = ', '.join('?' * 14)
                q = '({})'.format(q.replace('?', '%s'))

                key_str = get_key_str()
                insert_command = """INSERT INTO reaction
                ({0}) VALUES %s;""".format(key_str)

                execute_values(cur=cur, sql=insert_command,
                               argslist=reaction_values,
                               template=q, page_size=block_size)

                key_str = get_key_str('reaction_system')
                insert_command = """INSERT INTO reaction_system
                ({0}) VALUES %s ON CONFLICT DO NOTHING;""".format(key_str)

                execute_values(cur=cur, sql=insert_command,
                               argslist=reaction_system_values,
                               page_size=1000)
                con.commit()

                t2 = time.time()
                dt = t2 - t1
                t_av = (t_av * i + dt) / (i + 1)

                self.stdout.write(
                    '  Finnished Block {0} / {1} in {2} sec \n'
                    .format(block_id + 1, n_blocks, dt))
                self.stdout.write(
                    '    Completed transfer of {0} reactions. \n'
                    .format(Ncat - Ncat0))
                self.stdout.write('    Estimated time left: {0} sec \n'
                                  .format(t_av * (n_blocks - block_id - 1)))

            self.stdout.write('  Completed transfer of reactions\n')

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
              sites=None, reaction_energy=None):
        con = self.connection or self._connect()
        self._initialize(con)
        cur = con.cursor()
        keys = 'pub_id, chemical_composition,  reactants, products'
        values = [pub_id, chemical_composition, reactants, products]
        placeholder = """'{}', '{}', '{}', '{}'"""
        if sites:
            placeholder += ", '{}'"
            keys += ', sites'
            values.append(sites)
        if reaction_energy:
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


def get_key_list(table='reaction', start_index=0):
    key_list = {'reaction': ['id', 'chemical_composition',
                             'surface_composition',
                             'facet', 'sites', 'coverages', 'reactants',
                             'products', 'reaction_energy',
                             'activation_energy', 'dft_code',
                             'dft_functional', 'username', 'pub_id'],
                'publication': ['id', 'pub_id', 'title', 'authors', 'journal',
                                'volume', 'number', 'pages', 'year',
                                'publisher', 'doi', 'tags'],
                'reaction_system': ['name', 'energy_correction',
                                    'ase_id', 'id'],
                'publication_system': ['ase_id, pub_id'],
                'systems': ['id', 'unique_id', 'ctime', 'mtime', 'username',
                            'numbers', 'positions', 'cell', 'pbc',
                            'initial_magmoms', 'initial_charges', 'masses',
                            'tags', 'momenta', 'constraints', 'calculator',
                            'calculator_parameters', 'energy', 'free_energy',
                            'forces', 'stress', 'dipole', 'magmoms', 'magmom',
                            'charges', 'key_value_pairs', 'data', 'natoms',
                            'fmax', 'smax', 'volume', 'mass', 'charge']}

    return key_list[table][start_index:]


def get_key_str(table='reaction', start_index=0):
    key_str = """, """.join(get_key_list(table, start_index))

    return key_str


def get_value_list(values, start_index=0):
    value_list = []
    for v in values:
        if PY2:  # python 2
            if isinstance(v, unicode):
                v = v.encode('utf8', 'ignore')
        value_list += [v]
    return value_list[start_index:]


def get_value_str(values, start_index=0):
    values = get_value_list(values, start_index)
    value_str = "'{0}'".format(values[0])
    for v in values[1:]:
        if v is None or v == '':
            value_str += ", {0}".format('NULL')
        elif isinstance(v, str):
            value_str += ", '{0}'".format(v)
        else:
            value_str += ", {0}".format(v)
    return value_str
