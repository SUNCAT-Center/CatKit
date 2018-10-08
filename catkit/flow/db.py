from .. import Gratoms
from .. import gen
from . import utils
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import JSONB
import sqlalchemy as sqa
import passlib.apps
import operator
import warnings
import numpy as np
import json
import os
try:
    from ckutil.aflow import get_prototype_tag
except(ImportError):
    pass
Base = declarative_base()


class CalculationRequired(Exception):
    """Base class for unhandled calculation submission requirement."""
    pass


class Connect():
    """A class for accessing a temporary SQLite database. This
    function works as a context manager and should be used as follows:

    with Connect() as db:
        (Perform operation here)
    """

    def __init__(self, engine=None, user=None):
        """Initialize the database."""
        if engine is None:
            engine = os.environ.get('CATFLOW_DB')
        self.engine = sqa.create_engine(engine)
        self.session = sqa.orm.sessionmaker(bind=self.engine)
        self.meta = Base.metadata
        self.meta.create_all(self.engine)
        self.user = user

    def __enter__(self):
        """This function is automatically called whenever the class
        is used together with a 'with' statement.
        """
        self.cursor = self.session()

        return self

    def __exit__(self, type, value, tb):
        """Upon exiting the 'with' statement, __exit__ is called."""
        self.cursor.commit()
        self.cursor.close()

    def _reset(self, safety=False):
        """Drop the entire database schema.
        WARNING: This will remove all existing data from the database.

        Parameters
        ----------
        safety : bool
            A flag to prevent accidental deletion.
        """
        if safety:
            self.meta.drop_all(self.engine)
            print('Database reset')
        else:
            warnings.warn('You are about to reset the entire database! '
                          'This will delete all existing information. '
                          'If this is correct, set safety=True.')

    def structure_to_atoms(
            self,
            structure,
            parameters=None):
        """Return an atoms object for a given structure from the
        database.

        Parameters
        ----------
        structure : Structure object
            A structure from the CatFlow database.
        parameters : dict
            Parameters used to perform the calculation. If None, they
            will be retrieved from the database.

        Returns
        -------
        atoms : Gratoms object
            Atomic structure.
        """
        if parameters is None:
            calculator_name, parameters = self.cursor.query(
                Calculator.name, Calculator.parameters).\
                filter(Calculator.id == structure.calculator_id).one()
        else:
            calculator_name = parameters.pop('calculator_name')
        calculator = utils.str_to_class(calculator_name)

        atoms = Gratoms(
            numbers=structure.numbers,
            positions=structure.positions,
            pbc=structure.pbc,
            cell=structure.cell)

        results = {}
        for prop_name in utils.supported_properties:
            prop = getattr(structure, prop_name)
            if isinstance(prop, list):
                prop = np.array(prop)
            results[prop_name] = prop

        calculator = calculator(atoms, **parameters)
        calculator.set_results(results)

        return atoms

    def get_structure_trajectory(self, structure, toatoms=False):
        """Return the complete trajectory of a structure from the database.
        The trajectory is always returned in original relaxation step order.

        Parameters
        ----------
        structure : Structure object | int
            A structure from the CatFlow database or the table ID of the
            desired structure.
        toatoms : bool
            Whether to convert the list of structures to a list of atoms.

        Returns
        -------
        trajectory : list of (Structure objects | Atoms objects)
            The trajectory of calculations in relaxation step order.
        """
        check_initial = True
        if isinstance(structure, int):
            structure_id = structure
            check_initial = False
        else:
            structure_id = structure.id

        included_parts = self.cursor.query(Structure)\
                                    .filter(Structure.id == structure_id)\
                                    .cte(name='trajectroy', recursive=True)

        if check_initial and structure.parent_id is None:
            included_parts = included_parts.union_all(
                self.cursor.query(Structure)
                .filter(Structure.parent_id == included_parts.c.id))
            trajectory = self.cursor.query(included_parts).all()
        else:
            included_parts = included_parts.union_all(
                self.cursor.query(Structure)
                .filter(Structure.id == included_parts.c.parent_id))
            trajectory = self.cursor.query(included_parts)[::-1]

        if toatoms:
            images = []
            for structure in trajectory:
                images += [self.structure_to_atoms(structure)]
            trajectory = images

        return trajectory

    def get_matching_calculator(self, parameters, ignored=None):
        """Search for a calculator which matches a set of calculator
        parameters. This will raise an exception if more than one
        calculator is found.

        Parameters
        ----------
        parameters : dict
            Calculator input parameters to look for match.

        Returns
        -------
        calculator : Calculator object
            Catflow calculator which matches the requested parameters.
        """
        parameters = parameters.copy()
        name = parameters.pop('calculator_name')
        query = self.cursor.query(Calculator)\
                           .filter(Calculator.name == name)\
                           .filter(Calculator.parameters.contains(parameters))\
                           .filter(
                               Calculator.parameters.contained_by(parameters))
        calculator = query.one_or_none()

        return calculator

    def get_similar_calculators(self, parameters=None, comparator='ge'):
        """Search for a calculator which has convergence settings which
        are compared to the requested parameters by the requested
        comparator.

        Currently, only the decaf.Espresso calculator is supported.

        Parameters
        ----------
        parameters : dict
            Calculator input parameters to compare against.
        comparator : str
            A comparator type to use for parameter comparison:
            https://docs.python.org/2/library/operator.html

        Returns
        -------
        calculators : list of Calculator objects
            Catflow calculators with higher convergence criteria than then
            those provided.
        """
        parameters = parameters.copy()
        parameters.pop('calculator_name', None)
        convergence_criteria = [
            'kspacing', 'ecutwfc', 'conv_thr',
            'forc_conv_thr', 'press_conv_thr']

        criteria = {}
        for k in convergence_criteria:
            value = parameters.pop(k, None)
            if value:
                criteria[k] = value

        # Non-convergence criteria must match exactly
        query = self.cursor.query(Calculator)\
                           .filter(Calculator.name == 'decaf.Espresso')\
                           .filter(Calculator.parameters.contains(parameters))

        calculators = []
        for calculator in query:
            compare = getattr(operator, comparator)
            checks_passed = [
                compare(calculator.kspacing**-1,
                        criteria.get('kspacing', 1)**-1),
                compare(calculator.energy_cutoff,
                        criteria.get('ecutwfc', 1e12)),
                compare(calculator.scf_convergence,
                        criteria.get('conv_thr', 1e-4)),
                compare(calculator.force_convergence,
                        criteria.get('forc_conv_thr', 0.05)),
                compare(calculator.stress_convergence,
                        criteria.get('press_conv_thr', 0.5))]

            if np.all(checks_passed):
                calculators += [calculator]

        return calculators

    def request_bulk_entry(self, atoms, parameters, auto_submit=False):
        """Requests a bulk entry which matches a given atomic structures
        symmetry and set of calculator parameters. If a minimum structure
        is found, return it. Otherwise, prepare for a submission step.

        Parameters
        ----------
        atoms : Atoms object | str
            Initial guess of atomic structure to search for a symmetry match
            for. Will also search for matching prototype tag.
        parameters : dict
            Calculator parameters required to match for relaxed structure.
        auto_submit : bool
            Automatically perform a calculation submission if the requested
            structure is not preasent.

        Returns
        -------
        atoms : Gratoms object
            A matching relaxed bulk structure from the database.
        """
        # A.1) We need to know if this calculator exists.
        calculator = self.get_matching_calculator(parameters)

        if calculator is None:
            if not auto_submit or isinstance(atoms, str):
                raise CalculationRequired('No bulk structure match')
            # TODO: a more intellegent framework might use a different
            # calculator match as a good initial guess here.

            # B.1) Submit calculation.

        # A.2) We need a prototype tag to match against the database.
        if isinstance(atoms, str):
            tag = atoms
        else:
            tag, atoms = get_prototype_tag(atoms, tol=1e-2)
            atoms.info['prototype_tag'] = tag

        # A.3) Query for all similar structures with matching calculator
        query = self.cursor.query(Bulk.id, Bulk.prototype_tag, Structure).\
                    filter(Bulk.initial_prototype_tag == tag).\
                    join(Structure, Bulk.structure_id == Structure.id).\
                    filter(Structure.calculator_id == calculator.id)

        # A.4) If multiple structures exist, there could be multiple minima
        # User intervention will likely be required.
        # For now, lets assume there is one minima at most.
        previous_structure = query.one_or_none()
        if previous_structure is None:
            if not auto_submit:
                raise CalculationRequired('No bulk structure match')

            # B.1) Submit calculation.
            return
        elif previous_structure[1] is None:
            # A.5) No final structure. This calculation has already been submitted.
            # Print the current status.
            bulk_id, _, workflow, user = self.cursor.query(
                Bulk.id, Structure.worflow_id, Workflow, User).\
                filter(Bulk.id == previous_structure[0]).\
                join(Structure, Bulk.initial_structure_id == Structure.id).\
                join(Workflow).join(User).one()
            print(workflow)
            return workflow

        # A.6) If found, return the atoms object
        tag, structure = previous_structure
        atoms = self.structure_to_atoms(structure)

        return atoms

    def submit_bulk_entry(self, atoms, parameters):
        """Submit a workflow to an external workflow management software.
        Updates the CatFlow database with initial entries to indicate submission.

        Only integration with Fireworks is currently supported.

        Parameters
        ----------
        atoms : Atoms object
            The initial configuration of atoms to use for relaxation.
        parameters : dict
            The calculator parameters to use for the relaxation.
        """
        # Redundent for safety
        calculator = self.get_matching_calculator(parameters)

        if calculator is None:
            calculator = Calculator(parameters)
            self.cursor.add(calculator)

        user = self.cursor.query(User).filter(User.username = self.user)
        if user is None:
            raise RuntimeError('No workflow username match for {}. '.format(self.user)
                               'Please provide a user to Connect(). '
                               'New users must register with Connect().register_user()')

        # If the tag exists, assume the system is standardized
        tag = atoms.info.get('prototype_tag')
        if not tag:
            tag, atoms = get_prototype_tag(atoms, tol=1e-2)

        initial_structure = Structure(atoms, calculator=calculator)
        # Needed to resolve IDs
        self.cursor.commit()

        # Now Laminar

        # Need to submit the firework next and collect its ID
        workflow = Workflow(
            user, workflow_id, workflow_name='bulk_relaxation')

        initial_structure.update().values(workflow=workflow)

    def update_bulk_entry(self, trajectory, parameters, bulk_id, structure_id):
        """Add a bulk relaxation trajectory to the database. This requires
        the calculator parameters, the calculator being used, and the
        images from the trajectory.

        Parameters
        ----------
        trajectory : list of Atoms object
            Relaxation trajectory for a bulk calculation given in
            correct relaxation order.
        parameters : dict
            Calculator parameters required to match for relaxed structure.
        """
        try:
            atoms = self.request_bulk_entry(trajectory[0], parameters)
            warnings.warn('Matching calculation found, aborting')
            return
        except(CalculationRequired):
            pass

        calculator = self.get_matching_calculator(parameters)
        if calculator is None:
            calculator = Calculator(parameters)
            self.cursor.add(calculator)

        init_prototype_tag, standard_atoms = get_prototype_tag(
            trajectory[0], tol=1e-2)
        initial_structure = Structure(
            standard_atoms,
            calculator=calculator,
            results=trajectory[0].calc.results)
        self.cursor.add(initial_structure)

        for i, atoms in enumerate(trajectory[1:-1]):
            structure = Structure(
                atoms, calculator=calculator,
                parent=structure if i > 0 else initial_structure)
            self.cursor.add(structure)

        tag, standard_atoms = get_prototype_tag(trajectory[-1], tol=1e-2)
        structure = Structure(
            standard_atoms,
            calculator=calculator,
            parent=structure,
            results=trajectory[-1].calc.results)
        self.cursor.add(structure)

        if tag != init_prototype_tag:
            # TODO: Need to do a difference check here to ensure structure
            # is the same. For now, lets assume any prototype match is
            # for an identical minima.
            query = self.cursor.query(Bulk.prototype_tag, Structure).\
                        filter(Bulk.prototype_tag == tag).join(Structure).\
                        filter(Structure.calculator_id == calculator.id)
            previous_structure = query.one_or_none()

            if previous_structure is not None:
                # We have converged to a new prototype, need to
                # point to the previous relaxed structure.
                structure = previous_structure[-1]

        bulk_prototype = Bulk(
            init_prototype_tag, initial_structure, tag, structure=structure)

        self.cursor.add(bulk_prototype)

    def _add_bulk_entry(self, trajectory, parameters):
        """Add a bulk relaxation trajectory to the database. This requires
        the calculator parameters, the calculator being used, and the
        images from the trajectory.

        WARNING: Before being entered into the database, the structures will be
        standardized and the prototype tags recalculated. This is potentially
        quite dangerous as it assumes the structure is mostly unchanged. This
        is a very bad assumption, so this should only be used for testing.

        Parameters
        ----------
        trajectory : list of Atoms object
            Relaxation trajectory for a bulk calculation given in
            correct relaxation order.
        parameters : dict
            Calculator parameters required to match for relaxed structure.
        """
        try:
            atoms = self.request_bulk_entry(trajectory[0], parameters)
            warnings.warn('Matching calculation found, aborting')
            return
        except(CalculationRequired):
            pass

        calculator = self.get_matching_calculator(parameters)

        if calculator is None:
            calculator = Calculator(parameters)
            self.cursor.add(calculator)

        init_prototype_tag, standard_atoms = get_prototype_tag(
            trajectory[0], tol=1e-2)
        initial_structure = Structure(
            standard_atoms,
            calculator=calculator,
            results=trajectory[0].calc.results)
        self.cursor.add(initial_structure)

        for i, atoms in enumerate(trajectory[1:-1]):
            structure = Structure(
                atoms, calculator=calculator,
                parent=structure if i > 0 else initial_structure)
            self.cursor.add(structure)

        tag, standard_atoms = get_prototype_tag(trajectory[-1], tol=1e-2)
        structure = Structure(
            standard_atoms,
            calculator=calculator,
            parent=structure,
            results=trajectory[-1].calc.results)
        self.cursor.add(structure)

        if tag != init_prototype_tag:
            # TODO: Need to do a difference check here to ensure structure
            # is the same. For now, lets assume any prototype match is
            # for an identical minima.
            query = self.cursor.query(Bulk.prototype_tag, Structure).\
                        filter(Bulk.prototype_tag == tag).join(Structure).\
                        filter(Structure.calculator_id == calculator.id)
            previous_structure = query.one_or_none()

            if previous_structure is not None:
                # We have converged to a new prototype, need to
                # point to the previous relaxed structure.
                structure = previous_structure[-1]

        bulk_prototype = Bulk(
            init_prototype_tag, initial_structure, tag, structure=structure)

        self.cursor.add(bulk_prototype)


class User(Base):
    __tablename__ = 'user'
    id = sqa.Column(sqa.Integer, primary_key=True)

    first_name = sqa.Column(sqa.String, nullable=False)
    last_name = sqa.Column(sqa.String, nullable=False)

    host = sqa.Column(sqa.String, nullable=False)
    name = sqa.Column(sqa.String, nullable=False)
    username = sqa.Column(sqa.String, nullable=False)
    password_hash = sqa.Column(sqa.String, nullable=False)

    sqa.UniqueConstraint(host, name, username)

    def __init__(self, host, name, username, password
                 first_name, last_name):
        self.name = name
        self.host = host
        self.username = username
        self.password_hash = password_hash

        self.first_name = first_name
        self.last_name = last_name


class Workflow(Base):
    __tablename__ = 'workflow'
    id = sqa.Column(sqa.Integer, primary_key=True)

    user_id = sqa.Column(sqa.Integer, sqa.ForeignKey('user.id'))
    user = sqa.orm.relationship('User', uselist=False)
    workflow_id = sqa.Column(sqa.Integer, nullable=False)
    workflow_name = sqa.Column(sqa.String, nullable=False)

    sqa.UniqueConstraint(user_id, workflow_id)

    def __init__(self, user, worflow_name, workflow_id):
        self.user = user
        self.workflow_name = workflow_name
        self.workflow_id = workflow_id


class Calculator(Base):
    __tablename__ = 'calculator'
    id = sqa.Column(sqa.Integer, primary_key=True)

    name = sqa.Column(sqa.String, nullable=False)
    kspacing = sqa.Column(sqa.Float)
    energy_cutoff = sqa.Column(sqa.Integer)
    scf_convergence = sqa.Column(sqa.Float)
    force_convergence = sqa.Column(sqa.Float)
    stress_convergence = sqa.Column(sqa.Float)
    parameters = sqa.Column(JSONB)

    def __init__(self, parameters):
        parameters = parameters.copy()
        name =  parameters.pop('calculator_name', None)
        if not name:
            raise ValueError("parameters must contain 'calculator_name'")

        self.name = name

        # Special kwargs support for decaf-espresso
        if name == 'decaf.Espresso':
            # Energy convergence is not independent of
            # force convergence in decaf-espresso
            self.kspacing = parameters.get('kspacing', 1)
            self.energy_cutoff = parameters.get('ecutwfc', None)
            self.scf_convergence = parameters.get('conv_thr', 1e-4)
            self.force_convergence = parameters.get('forc_conv_thr', 0.05)
            self.stress_convergence = parameters.get('press_conv_thr', 0.5)

        self.parameters = parameters

    def __repr__(self):
        prompt = 'Calculator: ({}) {}\n'.format(self.id, self.name)
        prompt += json.dumps(self.parameters,
                             separators=(',', ': '),
                             indent=4, sort_keys=True)
        return prompt


class Structure(Base):
    __tablename__ = 'structure'
    id = sqa.Column(sqa.Integer, primary_key=True)

    numbers = sqa.Column(sqa.ARRAY(sqa.Integer), nullable=False)
    positions = sqa.Column(sqa.ARRAY(sqa.Float), nullable=False)
    cell = sqa.Column(sqa.ARRAY(sqa.Float), nullable=False)
    pbc = sqa.Column(sqa.ARRAY(sqa.Integer), nullable=False)

    # Optionals
    tags = sqa.Column(sqa.ARRAY(sqa.Integer))
    constraints = sqa.Column(sqa.ARRAY(sqa.Integer))

    # Results
    energy = sqa.Column(sqa.Float)
    forces = sqa.Column(sqa.ARRAY(sqa.Float))
    stress = sqa.Column(sqa.ARRAY(sqa.Float))
    charges = sqa.Column(sqa.ARRAY(sqa.Float))
    magmom = sqa.Column(sqa.Float)
    magmoms = sqa.Column(sqa.ARRAY(sqa.Float))

    calculator_id = sqa.Column(sqa.Integer, sqa.ForeignKey('calculator.id'))
    calculator = sqa.orm.relationship('Calculator', uselist=False)
    parent_id = sqa.Column(sqa.Integer, sqa.ForeignKey('structure.id'))
    parent = sqa.orm.relationship(
        'Structure', uselist=False, remote_side=[id],
        cascade="all, delete-orphan", single_parent=True)

    workflow_id = sqa.Column(sqa.Integer, sqa.ForeignKey('workflow.id'))
    workflow = sqa.orm.relationship('Workflow', uselist=False)

    def __init__(
            self,
            atoms,
            calculator=None,
            parent=None,
            workflow=None,
            results=None):
        self.numbers = atoms.numbers.tolist()
        self.positions = atoms.positions.tolist()
        self.cell = atoms.cell.tolist()
        self.pbc = atoms.pbc.astype(int).tolist()

        self.tags = None
        if any(atoms.get_tags()):
            self.tags = atoms.get_tags().tolist()

        # Read the constraints into an array
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

        # Read results
        if results is None and atoms.calc:
            results = atoms.calc.results
        for parameter in utils.supported_properties:
            result = results.get(parameter)
            if isinstance(result, np.ndarray):
                result = result.tolist()
            setattr(self, parameter, result)

        self.calculator = calculator
        self.parent = parent
        self.workflow = workflow


class Bulk(Base):
    __tablename__ = 'bulk'
    id = sqa.Column(sqa.Integer, primary_key=True)

    initial_prototype_tag = sqa.Column(sqa.String, nullable=False)
    initial_structure_id = sqa.Column(sqa.Integer, sqa.ForeignKey('structure.id'))
    initial_structure = sqa.orm.relationship(
        'Structure', uselist=False, foreign_keys=initial_structure_id)

    # Final state
    prototype_tag = sqa.Column(sqa.String)
    structure_id = sqa.Column(sqa.Integer, sqa.ForeignKey('structure.id'))
    structure = sqa.orm.relationship(
        'Structure', uselist=False, foreign_keys=structure_id)

    spacegroup = sqa.Column(sqa.Integer)
    wyckoff_positions = sqa.Column(sqa.String)
    composition = sqa.Column(sqa.String)

    def __init__(
            self,
            initial_prototype_tag,
            initial_structure=None,
            prototype_tag=None,
            structure=None):
        self.initial_prototype_tag = initial_prototype_tag
        self.initial_structure = initial_structure

        self.prototype_tag = prototype_tag
        if prototype_tag:
            details = prototype_tag.split('_')
            self.spacegroup = details[0]
            self.wyckoff_positions = details[1]
            self.composition = details[2]
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
