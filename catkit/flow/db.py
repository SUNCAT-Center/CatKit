from .. import Gratoms
from .. import gen
from . import Laminar
from . import utils
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import JSONB
import sqlalchemy as sqa
import numpy as np
import datetime
import operator
import warnings
import json
import ase
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

    def register_user(self, first_name, last_name,
                      host, name, username, password):
        """Register a users workflow credentials to the database.
        TODO: 2-way encryption for password protection.

        Currently, only Fireworks credentials are supported for automatic
        workflow submission.

        Parameters
        ----------
        first_name : str
            Users first name, for identification purposes.
        last_name : str
            Users last name, for identification purposes.
        host : str
            Server name of the fireworks database.
        name : str
            Table name of the fireowrks database.
        username : str
            Username of the fireworks database.
        password : str
            Password to access the fireworks database with.
        """
        user = User(first_name=first_name,
                    last_name=last_name,
                    host=host,
                    name=name,
                    username=username,
                    password=password)
        self.cursor.add(user)
        self.cursor.commit()

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
            self.submit_bulk_entry(atoms, parameters)
            return

        # A.2) We need a prototype tag to match against the database.
        if isinstance(atoms, str):
            tag = atoms
        else:
            tag, atoms = get_prototype_tag(atoms, tol=1e-2)
            atoms.info['prototype_tag'] = tag

        # A.3) Query for all similar structures with matching calculator
        query = self.cursor.query(Bulk.id, Bulk.prototype_tag, Structure).\
                    filter(Bulk.initial_prototype_tag == tag).\
                    join(Structure, Bulk.initial_structure_id == Structure.id).\
                    filter(Structure.calculator_id == calculator.id)

        # A.4) If multiple structures exist, there could be multiple minima
        # User intervention will likely be required.
        # For now, lets assume there is one minima at most.
        previous_structure = query.one_or_none()

        if previous_structure is None:
            if not auto_submit:
                raise CalculationRequired('No bulk structure match')

            # B.1) Submit calculation.
            self.submit_bulk_entry(atoms, parameters)
            return
        elif previous_structure[1] is None:
            # A.5) No final structure. This calculation has already been submitted.
            # Print the current status.
            bulk_id, _, workflow, user = self.cursor.query(
                Bulk.id, Structure.workflow_id, Workflow, User).\
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

        user = self.cursor.query(User).\
            filter(User.username == self.user).one_or_none()
        if user is None:
            raise RuntimeError(
                'No workflow username match for {}. Please provide a user '
                'to Connect(). New users must register with Connect().'
                'register_user()'.format(self.user))

        # If the tag exists, assume the system is standardized
        tag = atoms.info.get('prototype_tag')
        if not tag:
            tag, atoms = get_prototype_tag(atoms, tol=1e-2)

        initial_structure = Structure(atoms, calculator=calculator)
        self.cursor.add(initial_structure)

        bulk_prototype = Bulk(tag, initial_structure)
        self.cursor.add(bulk_prototype)
        self.cursor.commit()
        atoms.info['structure_id'] = initial_structure.id
        atoms.info['bulk_id'] = bulk_prototype.id

        # Now Laminar submission
        flow = Laminar(
            host=user.host,
            name=user.name,
            username=user.username,
            password=user.password_hash)
        workflow_id = flow.bulk_relaxation(atoms, parameters)

        # Need to submit the firework next and collect its ID
        workflow = Workflow(
            user, workflow_id, workflow_name='bulk_relaxation')

        initial_structure.workflow = workflow
        self.cursor.add(initial_structure)
        self.cursor.commit()

    def update_bulk_entry(self, trajectory, bulk_id=None, prototype_tag=None):
        """Update a bulk calcualtion which has been submitted with
        submit_bulk_entry(). This involves updating the information for an
        existing bulk entry.

        Parameters
        ----------
        trajectory : list of Atoms object
            Relaxation trajectory for a bulk calculation given in
            correct relaxation order.
        bulk_id : int
            Bulk entry to be updated.
        prototype_tag : str
            Prototype tag classification for the final structure.
        """
        if bulk_id is None:
            bulk_id = trajectory[0].info['bulk_id']

        bulk, initial_structure, workflow, calculator = self.cursor.\
            query(Bulk, Structure, Workflow, Calculator).\
            filter(Bulk.id == bulk_id).\
            join(Structure, Bulk.initial_structure_id == Structure.id).\
            join(Workflow).join(Calculator).one()

        for parameter in utils.supported_properties:
            result = trajectory[0].calc.results.get(parameter)
            if isinstance(result, np.ndarray):
                result = result.tolist()
            setattr(initial_structure, parameter, result)

        for i, atoms in enumerate(trajectory[1:-1]):
            structure = Structure(
                atoms,
                calculator=calculator,
                parent=structure if i > 0 else initial_structure,
                workflow=workflow)
            self.cursor.add(structure)

        structure = Structure(
            atoms,
            calculator=calculator,
            parent=structure,
            workflow=workflow)
        self.cursor.add(structure)

        if prototype_tag is None:
            prototype_tag = trajectory[-1].info['prototype_tag']

        bulk.prototype_tag = prototype_tag
        bulk.structure = structure

        details = prototype_tag.split('_')
        bulk.spacegroup = details[0]
        bulk.wyckoff_positions = details[1]
        bulk.composition = details[2]

        workflow.status = 'Completed'
        self.cursor.commit()

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

    def __init__(self, host, name, username, password,
                 first_name, last_name):
        self.name = name
        self.host = host
        self.username = username
        self.password_hash = password

        self.first_name = first_name
        self.last_name = last_name

    def __repr__(self):
        prompt = 'User: ({}) {} {}\n'.format(
            self.id, self.first_name, self.last_name)
        prompt += 'server: {}@{}/{}'.format(
            self.username, self.host, self.name)
        return prompt


class Workflow(Base):
    __tablename__ = 'workflow'
    id = sqa.Column(sqa.Integer, primary_key=True)

    user_id = sqa.Column(sqa.Integer, sqa.ForeignKey('user.id'))
    user = sqa.orm.relationship('User', uselist=False)
    workflow_id = sqa.Column(sqa.Integer, nullable=False)
    workflow_name = sqa.Column(sqa.String, nullable=False)
    status = sqa.Column(sqa.String, nullable=False)
    submission_date = sqa.Column(sqa.DateTime, nullable=False)

    sqa.UniqueConstraint(user_id, workflow_id)

    def __init__(self, user, workflow_id, workflow_name):
        self.user = user
        self.workflow_name = workflow_name
        self.workflow_id = workflow_id
        self.status = 'Submitted'
        self.submission_date = datetime.datetime.now().\
            strftime("%Y-%m-%d %H:%M:%S")

    def __repr__(self):
        prompt = 'Workflow: ({}) {}\n'.format(self.id, self.workflow_name)
        prompt += 'Status - {}'.format(self.status)
        if self.status == 'Submitted':
            prompt += '\n{}\n'.format(self.submission_date)
            prompt += '\n{}'.format(self.user)
        return prompt



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

        if results:
            for parameter in utils.supported_properties:
                result = results.get(parameter)
                if isinstance(result, np.ndarray):
                    result = result.tolist()
                setattr(self, parameter, result)

        self.calculator = calculator
        self.parent = parent
        self.workflow = workflow

    def __repr__(self):
        cid, pid, wid = 'None', 'None', 'None'
        if self.calculator:
            cid = self.calculator.id
        if self.parent:
            pid = self.parent.id
        if self.workflow:
            wid = self.workflow.id
        prompt = 'Structure : ({:<4})  Calculator: ({:<4})\n'.format(self.id, cid)
        prompt += 'Parent    : ({:<4})  Workflow  : ({:<4})\n\n'.format(pid, wid)

        cell = np.round(self.cell, 3)
        prompt += '| {:^20} | {:^10} \n |-\n'.format('Cell', 'pbc')
        for i in range(3):
            prompt += '| {:<8} {:<8} {:<8} | {} |\n'.format(*cell[i], self.pbc[i])

        prompt += '\n| Sym | {:^26} |\n'.format('Positions')
        prompt += '|-----|' + '-' * 28 + '|\n'

        symbols = np.array(ase.data.chemical_symbols)[self.numbers]
        positions = np.round(self.positions, 4)
        for i, n in enumerate(symbols):
            prompt += '|{:>4} |'.format(n)
            prompt += ' {:<8} {:<8} {:<8} |'.format(*positions[i])
        return prompt


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

    def __repr__(self):
        iid, fid = 'None', 'None'
        if self.initial_structure:
            iid = self.initial_structure.id
        if self.structure:
            fid = self.structure.id 
        prompt = 'Bulk: ({})\n'.format(self.id)
        prompt += 'Initial structure: ({:<4}) {}\n'.format(
            iid, self.initial_prototype_tag)
        prompt += 'Final structure:   ({:<4}) {}'.format(
            fid, self.prototype_tag)
        return prompt


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
