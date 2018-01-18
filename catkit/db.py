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

Base = declarative_base()


class Connect():
    """ A class for accessing a temporary SQLite database. This
    function works as a context manager and should be used as follows:

    with Connect() as db:
        (Perform operation here)

    """

    def __init__(
            self,
            engine='sqlite:///example.db'):
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

            structure_result = Structure_Result(
                calculator,
                structure,
                energy,
                stress)
            self.cursor.add(structure_result)

        for i, atom in enumerate(atoms):

            tmp_atom = Atom(
                structure,
                atom.number,
                atom.position,
                constraints[i],
                atom.magmom,
                atom.charge)
            self.cursor.add(tmp_atom)

            if atoms._calc and calculator:
                atom_result = Atom_Result(
                    calculator,
                    tmp_atom,
                    forces[i])
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

    atom = relationship(
        'Atom',
        uselist=False,
        back_populates='element')

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

    element = relationship(
        'Element',
        back_populates='atom')
    structure = relationship(
        'Structure',
        back_populates='atoms')
    atom_result = relationship(
        'Atom_Result',
        back_populates='atom')

    def __init__(
            self,
            structure,
            number,
            coordinates,
            constraints,
            magmom,
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

    atoms = relationship(
        'Atom',
        back_populates='structure')
    structure_result = relationship(
        'Structure_Result',
        back_populates='structure')
    entry = relationship(
        'Entry',
        back_populates='structure')

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

    entry = relationship(
        'Entry',
        back_populates='calculator')
    atom_result = relationship(
        'Atom_Result',
        uselist=False,
        back_populates='calculator')
    structure_result = relationship(
        'Structure_Result',
        uselist=False,
        back_populates='calculator')

    def __init__(
            self,
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

    calculator = relationship(
        'Calculator',
        back_populates='atom_result')
    atom = relationship(
        'Atom',
        back_populates='atom_result')

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

    calculator = relationship(
        'Calculator',
        back_populates='structure_result')
    structure = relationship(
        'Structure',
        back_populates='structure_result')

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

    structure = relationship(
        'Structure',
        back_populates='entry')
    calculator = relationship(
        'Calculator',
        back_populates='entry')
    trajectory = relationship(
        'Structure',
        secondary='trajectories')

    def __init__(
            self,
            structure,
            calculator,
            trajectory,
            natoms,
            energy,
            fmax,
            smax=None,
            search_keys=None
    ):
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
