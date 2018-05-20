from ase import data
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, ForeignKey
from sqlalchemy import Integer, String, Numeric, Boolean
from sqlalchemy.orm import relationship
from sqlalchemy.dialects.postgresql import JSONB
import numpy as np
import json
import sqlite3
from sqlite3 import IntegrityError
Base = declarative_base()


class Connect():
    """ A class for accessing a temporary SQLite database. This
    function works as a context manager and should be used as follows:

    with Connect() as db:
        (Perform operation here)
    """

    def __init__(self, engine='sqlite:///example.db'):
        """ The __init__ function is automatically called when the
        class is referenced.
        """

        self.engine = create_engine(engine)
        self.session = sessionmaker(bind=self.engine)

    def __enter__(self):
        """ This function is automatically called whenever the class
        is used together with a 'with' statement.
        """

        Base.metadata.create_all(self.engine)
        self.cursor = self.session()

        return self

    def __exit__(self, type, value, tb):
        """ Upon exiting the 'with' statement, __exit__ is called.
        """

        self.cursor.commit()
        self.cursor.close()

    def add_entry(self, images, search_keys=None):

        if isinstance(images, list):
            atoms = images[-1]
        else:
            atoms = images
            images = [atoms]

        if atoms._calc:

            parameters = atoms.info

            calculator = Calculator(
                name='QE',
                xc=parameters['xc'],
                kpoints=parameters['kpts'],
                energy_cutoff=parameters['pw'],
                parameters=parameters)

            self.cursor.add(calculator)

        trajectory = []
        for atoms in images:
            structure = self._add_structure(atoms, calculator)
            trajectory += [structure]

        entry = Entry(
            structure=trajectory[-1],
            calculator=calculator,
            trajectory=trajectory,
            natoms=len(atoms),
            energy=atoms.get_potential_energy(),
            fmax=abs(atoms.get_forces()).max(),
            smax=abs(atoms.get_stress()).max(),
            search_keys=json.dumps(search_keys),
        )

        self.cursor.add(entry)

    def _initialize_elements(self):
        for i in range(96):
            el = Element((i + 1))
            self.cursor.add(el)

    def _add_structure(self, atoms, calculator=None):

        cell = atoms.get_cell()
        pbc = atoms.get_pbc()
        structure = Structure(cell, pbc)
        self.cursor.add(structure)

        constraints = self._parse_constraint(atoms)

        if atoms._calc and calculator:

            energy = atoms.get_potential_energy()
            forces = atoms.get_forces()
            atoms.constraints = None
            stress = atoms.get_stress()

            structure_result = Structure_Result(calculator, structure, energy,
                                                stress)
            self.cursor.add(structure_result)

        for i, atom in enumerate(atoms):

            tmp_atom = Atom(structure, atom.number, atom.position,
                            constraints[i], atom.magmom, atom.charge)
            self.cursor.add(tmp_atom)

            if atoms._calc and calculator:
                atom_result = Atom_Result(calculator, tmp_atom, forces[i])
                self.cursor.add(atom_result)

        return structure

    def _parse_constraint(self, atoms):

        constraints = np.zeros((len(atoms), 3), dtype=int)
        for entry in atoms.constraints:
            d = entry.todict()
            name = d['name']

            if name == 'FixAtoms':
                i = d['kwargs']['indices']
                constraints[i] += 1

            elif name == 'FixCartesian':
                i = d['kwargs']['a']
                constraints[i] = d['kwargs']['mask']

        # Casting as string breaks for SQLite
        return constraints.astype(str)


class Element(Base):
    __tablename__ = 'element'

    id = Column(Integer, primary_key=True)

    symbol = Column(String(2), nullable=False)
    mass = Column(Numeric, nullable=False)
    covalent_radii = Column(Numeric, nullable=False)

    atom = relationship('Atom', uselist=False, back_populates='element')

    def __init__(self, i):
        self.id = int(i)
        self.symbol = data.chemical_symbols[i]
        self.mass = data.atomic_masses[i]
        self.covalent_radii = data.covalent_radii[i]


class Atom(Base):
    __tablename__ = 'atom'

    id = Column(Integer, primary_key=True)

    element_id = Column(Integer, ForeignKey('element.id'))
    structure_id = Column(Integer, ForeignKey('structure.id'))
    x_coordinate = Column(Numeric, nullable=False)
    y_coordinate = Column(Numeric, nullable=False)
    z_coordinate = Column(Numeric, nullable=False)
    x_fixed = Column(Boolean)
    y_fixed = Column(Boolean)
    z_fixed = Column(Boolean)
    magmom = Column(Numeric)
    charge = Column(Numeric)

    element = relationship('Element', back_populates='atom')
    structure = relationship('Structure', back_populates='atoms')
    atom_result = relationship('Atom_Result', back_populates='atom')

    def __init__(self, structure, number, coordinates, constraints, magmom,
                 charge):
        self.structure = structure
        self.element_id = int(number)
        self.x_coordinate = coordinates[0]
        self.y_coordinate = coordinates[1]
        self.z_coordinate = coordinates[2]
        self.x_fixed = constraints[0]
        self.y_fixed = constraints[1]
        self.z_fixed = constraints[2]
        self.magmom = magmom
        self.charge = charge


class Structure(Base):
    __tablename__ = 'structure'

    id = Column(Integer, primary_key=True)

    x1_cell = Column(Numeric, nullable=False)
    x2_cell = Column(Numeric, nullable=False)
    x3_cell = Column(Numeric, nullable=False)
    y1_cell = Column(Numeric, nullable=False)
    y2_cell = Column(Numeric, nullable=False)
    y3_cell = Column(Numeric, nullable=False)
    z1_cell = Column(Numeric, nullable=False)
    z2_cell = Column(Numeric, nullable=False)
    z3_cell = Column(Numeric, nullable=False)
    pbc = Column(Integer, nullable=False)

    atoms = relationship('Atom', back_populates='structure')
    structure_result = relationship(
        'Structure_Result', back_populates='structure')
    entry = relationship('Entry', back_populates='structure')

    def __init__(self, cell, pbc):
        self.x1_cell = cell[0][0]
        self.x2_cell = cell[0][1]
        self.x3_cell = cell[0][2]
        self.y1_cell = cell[1][0]
        self.y2_cell = cell[1][1]
        self.y3_cell = cell[1][2]
        self.z1_cell = cell[2][0]
        self.z2_cell = cell[2][1]
        self.z3_cell = cell[2][2]
        self.pbc = int(sum(pbc))


class Calculator(Base):
    __tablename__ = 'calculator'

    id = Column(Integer, primary_key=True)

    name = Column(String)
    xc = Column(String)
    x_kpoints = Column(Integer)
    y_kpoints = Column(Integer)
    z_kpoints = Column(Integer)
    energy_cutoff = Column(Numeric)
    parameters = Column(JSONB)

    entry = relationship('Entry', back_populates='calculator')
    atom_result = relationship(
        'Atom_Result', uselist=False, back_populates='calculator')
    structure_result = relationship(
        'Structure_Result', uselist=False, back_populates='calculator')

    def __init__(self,
                 name,
                 xc=None,
                 kpoints=None,
                 energy_cutoff=None,
                 parameters=None):
        self.name = name
        self.xc = xc
        self.energy_cutoff = energy_cutoff
        self.parameters = parameters

        if isinstance(kpoints, list):
            self.x_kpoints = kpoints[0]
            self.y_kpoints = kpoints[1]
            self.z_kpoints = kpoints[2]


class Atom_Result(Base):
    __tablename__ = 'atom_result'

    id = Column(Integer, primary_key=True)

    calculator_id = Column(Integer, ForeignKey('calculator.id'))
    atoms_id = Column(Integer, ForeignKey('atom.id'))
    x_force = Column(Numeric)
    y_force = Column(Numeric)
    z_force = Column(Numeric)

    calculator = relationship('Calculator', back_populates='atom_result')
    atom = relationship('Atom', back_populates='atom_result')

    def __init__(self, calculator, atom, forces=None):
        self.calculator = calculator
        self.atom = atom

        if forces is None:
            return

        self.x_force = forces[0]
        self.y_force = forces[1]
        self.z_force = forces[2]


class Structure_Result(Base):
    __tablename__ = 'structure_result'

    id = Column(Integer, primary_key=True)

    calcualtor_id = Column(Integer, ForeignKey('calculator.id'))
    structure_id = Column(Integer, ForeignKey('structure.id'))
    energy = Column(Numeric)
    xx_stress = Column(Numeric)
    yy_stress = Column(Numeric)
    zz_stress = Column(Numeric)
    yz_stress = Column(Numeric)
    xz_stress = Column(Numeric)
    xy_stress = Column(Numeric)

    calculator = relationship('Calculator', back_populates='structure_result')
    structure = relationship('Structure', back_populates='structure_result')

    def __init__(self, calculator, structure, energy, stress=None):
        self.calculator = calculator
        self.structure = structure
        self.energy = energy

        if stress is None:
            return

        self.xx_stress = stress[0]
        self.yy_stress = stress[1]
        self.zz_stress = stress[2]
        self.yz_stress = stress[3]
        self.xz_stress = stress[4]
        self.xy_stress = stress[5]


class Entry(Base):
    __tablename__ = 'entry'

    id = Column(Integer, primary_key=True)

    structure_id = Column(Integer, ForeignKey('structure.id'))
    calculator_id = Column(Integer, ForeignKey('calculator.id'))
    natoms = Column(Integer)
    energy = Column(Numeric)
    fmax = Column(Numeric)
    smax = Column(Numeric)
    search_keys = Column(JSONB)

    structure = relationship('Structure', back_populates='entry')
    calculator = relationship('Calculator', back_populates='entry')
    trajectory = relationship('Structure', secondary='trajectories')

    def __init__(self,
                 structure,
                 calculator,
                 trajectory,
                 natoms,
                 energy,
                 fmax,
                 smax=None,
                 search_keys=None):
        self.structure = structure
        self.calculator = calculator
        self.trajectory = trajectory
        self.natoms = natoms
        self.energy = energy
        self.fmax = fmax
        self.smax = smax
        self.search_keys = search_keys


class Trajectories(Base):
    __tablename__ = 'trajectories'

    id = Column(Integer, primary_key=True)
    structure_id = Column(Integer, ForeignKey('structure.id'))
    entry_id = Column(Integer, ForeignKey('entry.id'))


class FingerprintDB():
    """ A class for accessing a temporary SQLite database. This
    function works as a context manager and should be used as follows:

    with FingerprintDB() as fpdb:
        (Perform operation here)

    This syntax will automatically construct the temporary database,
    or access an existing one. Upon exiting the indentation, the
    changes to the database will be automatically committed.
    """

    def __init__(self, db_name='fingerprints.db', verbose=False):
        """ The __init__ function is automatically called when the
        class is referenced.

        Args:
            db_name (str): Name of the database fileto access. Will
            connect to 'fingerprints.db' by default.
            verbose (bool): Will print additional information.
        """
        self.db_name = db_name
        self.verbose = verbose

    def __enter__(self):
        """ This function is automatically called whenever the class
        is used together with a 'with' statement.
        """
        self.con = sqlite3.connect(self.db_name)
        self.c = self.con.cursor()
        self.create_table()

        return self

    def __exit__(self, type, value, tb):
        """ Upon exiting the 'with' statement, __exit__ is called.
        """
        self.con.commit()
        self.con.close()

    def create_table(self):
        """ Creates the database table framework used in SQLite.
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
        """ Enters a single ase-db image into the fingerprint database.

        This table can be expanded to contain atoms objects in the future.

        Args:
            d (object): An ase-db object which can be parsed.
            identity (str): An identifier of the users choice.

        Returns:
            int: The ase ID colleted for the ase-db object.
        """
        # ase-db ID with identity must be unique. If not, it will be skipped.
        try:
            self.c.execute("""INSERT INTO images (ase_id, identity)
            VALUES(?, ?)""", (d.id, identity))
        except (IntegrityError):
            if self.verbose:
                print('ASE ID with identifier already defined: {} {}'.format(
                    d.id, identity))

        return d.id

    def parameter_entry(self, symbol=None, description=None):
        """ A function for entering unique parameters into the database.

        Parameters
        ----------
            symbol (str): A unique symbol the entry can be referenced
            by. If None, the symbol will be the ID of the parameter
            as a string.
            description (str): A description of the parameter.
        """
        # If no symbol is provided, use the parameter ID (str).
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
        """ Get an array of integer values which correspond to the
        parameter IDs for a set of provided symbols. If no selection
        is provided, return all symbols.

        Parameters
        ----------
            selection (list): List of symbols in parameters
            table to be selected.
            display (bool): If True, print parameter descriptions.

        Returns:
            1-d array: Return the integer values of selected parameters.
        """
        # If no selection is made, return all parameters.
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

        return np.array(res).T[0].astype(int)

    def fingerprint_entry(self, ase_id, param_id, value):
        """ Enters a fingerprint value to the database for a
        given ase and parameter ID.

        Parameters
        ----------
            ase_id (int): The ase unique ID associated with an atoms object
            in the database.
            param_id (int or str): The parameter ID or symbol associated
            with and entry in the paramters table.
            value (float): The value of the parameter for the atoms object.
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
        """ Get the array of values associated with the provided parameters
        for each ase_id provided.

        Parameters
        ----------
            ase_id (list): The ase ID(s) associated with an atoms object in
            the database.
            params (list): List of symbols or int in parameters table to be
            selected.

        Returns
        -------
            n-d array: An array of values associated with the given
            parameters (a fingerprint) for each ase_id.
        """
        if isinstance(params, np.ndarray):
            params = params.tolist()

        if not params or isinstance(params[0], str):
            params = self.get_parameters(selection=params)
            psel = ','.join(params.astype(str))
        elif isinstance(params[0], int):
            psel = ','.join(np.array(params).astype(str))

        if ase_ids is None:
            cmd = """SELECT GROUP_CONCAT(value) FROM fingerprints
            JOIN images on fingerprints.image_id = images.iid
            WHERE param_id IN ({})
            GROUP BY ase_id
            ORDER BY images.iid""".format(psel)

        else:
            asel = ','.join(np.array(ase_ids).astype(str))

            cmd = """SELECT GROUP_CONCAT(value) FROM fingerprints
            JOIN images on fingerprints.image_id = images.iid
            WHERE param_id IN ({}) AND ase_id IN ({})
            GROUP BY ase_id""".format(psel, asel)

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        fingerprint = np.zeros((len(fetch), len(params)))
        for i, f in enumerate(fetch):
            fingerprint[i] = f[0].split(',')

        return fingerprint
