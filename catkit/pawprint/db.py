import numpy as np
import sqlite3
from sqlite3 import IntegrityError


class FingerprintDB():
    """A class for accessing a temporary SQLite database. This
    function works as a context manager and should be used as follows:

    with FingerprintDB() as fpdb:
        (Perform operation here)

    This syntax will automatically construct the temporary database,
    or access an existing one. Upon exiting the indentation, the
    changes to the database will be automatically committed.
    """

    def __init__(self, db_name='fingerprints.db', verbose=False):
        """Initialize the database.

        Parameters
        ----------
        db_name : str
            Name of the database file to access.
        verbose : bool
            Print additional information
        """
        self.db_name = db_name
        self.verbose = verbose

    def __enter__(self):
        """This function is automatically called whenever the class
        is used together with a 'with' statement.
        """
        self.con = sqlite3.connect(self.db_name)
        self.c = self.con.cursor()
        self.create_table()

        return self

    def __exit__(self, type, value, tb):
        """Upon exiting the 'with' statement, __exit__ is called."""
        self.con.commit()
        self.con.close()

    def create_table(self):
        """Creates the database table framework used in SQLite.
        This includes 3 tables: images, parameters, and fingerprints.

        The images table currently stores ase_id information and
        a unqiue string. This can be adapted in the future to support
        atoms objects.

        The parameters table stores a symbol (10 character maximum)
        for convenient reference and a description of the parameter.

        The fingerprints table holds a unique image and parmeter ID
        along with a float value for each. The ID pair must be unique.
        """
        self.c.execute("""CREATE TABLE IF NOT EXISTS images(
        iid INTEGER PRIMARY KEY AUTOINCREMENT,
        ase_id CHAR(32) UNIQUE NOT NULL,
        identity TEXT
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS parameters(
        pid INTEGER PRIMARY KEY AUTOINCREMENT,
        symbol CHAR(10) UNIQUE NOT NULL,
        description TEXT
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS fingerprints(
        entry_id INTEGER PRIMARY KEY AUTOINCREMENT,
        image_id INT NOT NULL,
        param_id INT NOT NULL,
        value REAL,
        FOREIGN KEY(image_id) REFERENCES images(image_id)
        ON DELETE CASCADE ON UPDATE CASCADE,
        FOREIGN KEY(param_id) REFERENCES parameters(param_id)
        ON DELETE CASCADE ON UPDATE CASCADE,
        UNIQUE(image_id, param_id)
        )""")

    def image_entry(self, d, identity=None):
        """Enters a single ase-db image into the fingerprint database.
        The ase-db ID with identity must be unique. If not, it will be skipped.

        This table can be expanded to contain atoms objects in the future.

        Parameters
        ----------
        d : ase-db object
            Database entry to parse.
        identity : str
            An identifier of the users choice.

        Returns
        -------
        ase_id : int
            The ase id collected.
        """
        ase_id = d.id

        try:
            self.c.execute("""INSERT INTO images (ase_id, identity)
            VALUES(?, ?)""", (ase_id, identity))
        except (IntegrityError):
            if self.verbose:
                print('ASE ID with identifier already defined: {} {}'.format(
                    ase_id, identity))

        return ase_id

    def parameter_entry(self, symbol=None, description=None):
        """Enters a unique parameter into the database.

        Parameters
        ----------
        symbol : str
            A unique symbol the entry can be referenced by. If None,
            the symbol will be the ID of the parameter
            as a string.
        description : str
            A description of the parameter.
        """
        if not symbol:
            self.c.execute("""SELECT MAX(pid) FROM parameters""")
            symbol = str(int(self.c.fetchone()[0]) + 1)

        # The symbol must be unique. If not, it will be skipped.
        try:
            self.c.execute("""INSERT INTO parameters (symbol, description)
            VALUES(?, ?)""", (symbol, description))
        except (IntegrityError):
            if self.verbose:
                print('Symbol already defined: {}'.format(symbol))

        # Each instance needs to be commited to ensure no overwriting.
        # This could potentially result in slowdown.
        self.con.commit()

    def get_parameters(self, selection=None, display=False):
        """Get an array of integer values which correspond to the
        parameter IDs for a set of provided symbols.

        Parameters
        ----------
        selection : list
            Symbols in parameters table to be selected. If no selection
            is made, return all parameters.
        display : bool
            Print parameter descriptions.

        Returns
        -------
        parameter_ids : array (n,)
            Integer values of selected parameters.
        """
        if not selection:
            self.c.execute("""SELECT pid, symbol, description
            FROM parameters""")
            res = self.c.fetchall()
        else:
            res = []
            for i, s in enumerate(selection):
                self.c.execute("""SELECT pid, symbol, description
                FROM parameters WHERE symbol = '{}'""".format(s))
                res += [self.c.fetchone()]

        if display:
            print('[ID ]: key    - Description')
            print('---------------------------')
            for r in res:
                print('[{0:^3}]: {1:<10} - {2}'.format(*r))

        parameter_ids = np.array(res).T[0].astype(int)

        return parameter_ids

    def fingerprint_entry(self, ase_id, param_id, value):
        """Enters a fingerprint value to the database for a given ase and
        parameter id.

        Parameters
        ----------
        ase_id : int
            The unique id associated with an atoms object in the database.
        param_id : int or str
            The parameter ID or symbol associated with and entry in the
            parameters table.
        value : float
            The value of the parameter for the atoms object.
        """
        # If parameter symbol is given, get the ID
        if isinstance(param_id, str):
            self.c.execute("""SELECT pid FROM parameters
            WHERE symbol = '{}'""".format(param_id))
            param_id = self.c.fetchone()

            if param_id:
                param_id = param_id[0]
            else:
                raise (KeyError, 'parameter symbol not found')

        self.c.execute("""SELECT iid FROM images
        WHERE ase_id = {}""".format(ase_id))
        image_id = self.c.fetchone()[0]

        self.c.execute("""INSERT INTO fingerprints (image_id, param_id, value)
        VALUES(?, ?, ?)""", (int(image_id), int(param_id), float(value)))

    def get_fingerprints(self, ase_ids=None, params=[]):
        """Get the array of values associated with the provided parameters
        for each ase_id.

        Parameters
        ----------
        ase_id : list
            The ase-id associated with an atoms object in the database.
        params : list
            List of symbols or int in parameters table to be selected.

        Returns
        -------
        fingerprint : array (n,)
            An array of values associated with the given parameters
            (a fingerprint) for each ase_id.
        """
        if isinstance(params, np.ndarray):
            params = params.tolist()

        if not params or isinstance(params[0], str):
            params = self.get_parameters(selection=params)
            psel = ','.join(params.astype(str))
        elif isinstance(params[0], int):
            psel = ','.join(np.array(params).astype(str))

        if ase_ids is None:
            cmd = """SELECT GROUP_CONCAT(IFNULL(value, 'nan')) FROM
            fingerprints JOIN images on fingerprints.image_id = images.iid
            WHERE param_id IN ({})
            GROUP BY ase_id
            ORDER BY images.iid""".format(psel)

        else:
            asel = ','.join(np.array(ase_ids).astype(str))

            cmd = """SELECT GROUP_CONCAT(IFNULL(value, 'nan')) FROM
            fingerprints JOIN images on fingerprints.image_id = images.iid
            WHERE param_id IN ({}) AND ase_id IN ({})
            GROUP BY ase_id""".format(psel, asel)

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        fingerprint = np.zeros((len(fetch), len(params)))
        for i, f in enumerate(fetch):
            fingerprint[i] = f[0].split(',')

        return fingerprint
