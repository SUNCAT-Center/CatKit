from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import JSONB
import sqlalchemy as sqa
import numpy as np
import json
Base = declarative_base()


class Calculator(Base):
    __tablename__ = 'calculator'
    id = sqa.Column(sqa.Integer, primary_key=True)

    name = sqa.Column(sqa.String)
    kpoints = sqa.Column(sqa.ARRAY(sqa.Integer), nullable=False)
    parameters = sqa.Column(JSONB)

    def __init__(self, name, parameters=None):
        self.name = name
        self.parameters = parameters
        if isinstance(parameters, dict) and 'kpts' in parameters:
            self.kpoints = parameters.get('kpts')
            del parameters['kpts']
            if isinstance(self.kpoints, np.ndarray):
                self.kpoints = self.kpoints.tolist()
            self.parameters = json.dumps(parameters)
        else:
            self.kpoints = [1, 1, 1]


class Structure(Base):
    __tablename__ = 'structure'
    id = sqa.Column(sqa.Integer, primary_key=True)

    numbers = sqa.Column(sqa.ARRAY(sqa.Integer), nullable=False)
    positions = sqa.Column(sqa.ARRAY(sqa.Numeric), nullable=False)
    cell = sqa.Column(sqa.ARRAY(sqa.Numeric), nullable=False)
    pbc = sqa.Column(sqa.ARRAY(sqa.Integer), nullable=False)

    # Optionals
    tags = sqa.Column(sqa.ARRAY(sqa.Integer))
    constraints = sqa.Column(sqa.ARRAY(sqa.Integer))

    # Results
    energy = sqa.Column(sqa.Numeric)
    forces = sqa.Column(sqa.ARRAY(sqa.Numeric))
    stress = sqa.Column(sqa.ARRAY(sqa.Numeric))
    charges = sqa.Column(sqa.ARRAY(sqa.Numeric))
    magmom = sqa.Column(sqa.Numeric)
    magmoms = sqa.Column(sqa.ARRAY(sqa.Numeric))

    calculator_id = sqa.Column(sqa.Integer, sqa.ForeignKey('calculator.id'))
    calculator = sqa.orm.relationship('Calculator')
    parent_id = sqa.Column(sqa.Integer, sqa.ForeignKey('structure.id'))
    parent = sqa.orm.relationship(
        'Structure', uselist=False, remote_side=[id],
        cascade="all, delete-orphan", single_parent=True)

    def __init__(self, atoms, calculator=None, parent=None):
        self.numbers = atoms.numbers.tolist()
        self.positions = atoms.positions.tolist()
        self.cell = atoms.cell.tolist()
        self.pbc = atoms.pbc.astype(int).tolist()

        self.tags = None
        if any(atoms.get_tags()):
            self.tags = atoms.get_tags().tolist()

        self.constraints = None
        if atoms.constraints:
            self.constraints = np.ones((len(atoms), 3), dtype=int)
            for cons in atoms.constraints:
                c = cons.todict()
                if c['name'] == 'FixAtoms':
                    self.constraints[c['kwargs']['indices']] = 0
                elif c['name'] == 'FixCartesian':
                    self.constraints[c['kwargs']['a']] = ~c['kwargs']['mask']
            self.constraints = self.constraints.tolist()

        parameters = [
            'energy', 'forces', 'stress', 'charges', 'magmom', 'magmoms']

        if atoms._calc:
            for parameter in parameters:
                result = atoms._calc.results.get(parameter)
                if isinstance(result, np.ndarray):
                    result = result.tolist()
                setattr(self, parameter, result)

        self.calculator = calculator
        self.parent = parent


class Bulk(Base):
    __tablename__ = 'bulk'
    id = sqa.Column(sqa.Integer, primary_key=True)

    composition = sqa.Column(sqa.String, nullable=False)
    spacegroup = sqa.Column(sqa.Integer, nullable=False)
    wyckoff_positions = sqa.Column(sqa.String, nullable=False)

    structure_id = sqa.Column(sqa.Integer, sqa.ForeignKey('structure.id'))
    structure = sqa.orm.relationship('Structure', uselist=False)

    # There can be multiple structures based on different calculators
    sqa.UniqueConstraint(
        composition, spacegroup, wyckoff_positions, structure_id)

    def __init__(
            self,
            composition,
            spacegroup,
            wyckoff_positions,
            structure=None):
        self.composition = composition
        self.spacegroup = spacegroup
        self.wyckoff_positions = wyckoff_positions
        self.structure = structure


class Surface(Base):
    __tablename__ = 'surface'
    id = sqa.Column(sqa.Integer, primary_key=True)

    connectivity = sqa.Column(sqa.ARRAY(sqa.Integer), nullable=False)
    surface_atoms = sqa.Column(sqa.ARRAY(sqa.Integer), nullable=False)
    bulk_id = sqa.Column(
        sqa.Integer, sqa.ForeignKey('bulk.id'), nullable=False)
    bulk = sqa.orm.relationship(
        'Bulk', uselist=False,
        cascade="all, delete-orphan", single_parent=True)

    structure_id = sqa.Column(sqa.Integer, sqa.ForeignKey('structure.id'))
    structure = sqa.orm.relationship('Structure', uselist=False)

    def __init__(
            self,
            bulk,
            connectivity,
            surface_atoms,
            structure=None):
        self.bulk = bulk
        self.connectivity = connectivity
        self.surface_atoms = surface_atoms
        self.structure = structure


class AdsorptionMolecules(Base):
    __tablename__ = 'adsorptionmolecules'
    adsorption_id = sqa.Column(
        sqa.Integer, sqa.ForeignKey('adsorption.id'), primary_key=True)
    molecule_id = sqa.Column(
        sqa.Integer, sqa.ForeignKey('molecule.id'), primary_key=True)
    adsorption_site = sqa.Column(sqa.Integer)
    orientation = sqa.Column(sqa.Integer)

    molecule = sqa.orm.relationship('Molecule')


class Molecule(Base):
    __tablename__ = 'molecule'
    id = sqa.Column(sqa.Integer, primary_key=True)

    connectivity = sqa.Column(sqa.ARRAY(sqa.Integer), nullable=False)

    structure_id = sqa.Column(sqa.Integer, sqa.ForeignKey('structure.id'))
    structure = sqa.orm.relationship('Structure', uselist=False)

    def __init__(self, connectivity, structure=None):
        if isinstance(connectivity, np.ndarray):
            connectivity = connectivity.tolist()
        self.connectivity = connectivity
        self.structure = structure


class Adsorption(Base):
    __tablename__ = 'adsorption'
    id = sqa.Column(sqa.Integer, primary_key=True)

    connectivity = sqa.Column(sqa.ARRAY(sqa.Integer))
    surface_id = sqa.Column(
        sqa.Integer, sqa.ForeignKey('surface.id'), nullable=False)
    surface = sqa.orm.relationship(
        'Surface', uselist=False,
        cascade="all, delete-orphan", single_parent=True)

    molecules = sqa.orm.relationship('AdsorptionMolecules')

    structure_id = sqa.Column(sqa.Integer, sqa.ForeignKey('structure.id'))
    structure = sqa.orm.relationship('Structure', uselist=False)

    def __init__(
            self,
            surface,
            molecules,
            connectivity,
            structure=None):
        self.surface = surface
        self.molecules = molecules
        self.connectivity = connectivity
        self.structure = structure
