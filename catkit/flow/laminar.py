from . import fwio
import fireworks
from ase.db.row import AtomsRow
from ase.db import connect
from netrc import netrc
import numpy as np


class Laminar():
    """Simple submission script helper for CatFlow.
    """

    def __init__(
            self,
            host,
            username=None,
            name=None,
            password=None):
        """Initialize a fireworks instance."""
        if username is None or name is None or password is None:
            username, name, password = netrc().authenticators(host)
        if username is None:
            raise ValueError('username, name, and password required.')

        launchpad = fireworks.LaunchPad(
            host=host,
            name=name,
            username=username,
            password=password)

        self.launchpad = launchpad

    def submit_relaxation(
            self,
            image,
            workflow_name,
            parameters=None,
            spec=None):
        """Run a relaxation of a given DB entry or atoms object. If a
        database object is used, the calculation will automatically store
        the keys and data for later retrieval.

        The entries uuid will also be stored and `data.calculator_parameters`
        will be used as the calculation parameters.

        Parameters
        ----------
        images : Atoms object | AtomsRow object
            ASE database entry or atoms object to relax.
        workflow_name : str
            Name of the fireworks calculation to be used.
        parameters : dict
            Calculation parameters to use. Will be pulled from
            a database entry `data.calculator_parameters`.
        spec : dict
            Additional fireworks specifications to pass to the database.
        """
        keys, data = {}, {}
        if isinstance(image, AtomsRow):
            atoms = image.toatoms()
            keys.update(image.key_value_pairs)
            keys.update({'uuid': image.unique_id})
            data.update(image.data)
        else:
            atoms = image

        calculator = parameters.pop('calculator_name', None)
        if calculator is None:
            raise ValueError("'calculator_name' missing from parameters.")

        if parameters is not None:
            atoms.info['calculator_parameters'] = parameters
        elif data.get('calculator_parameters'):
            atoms.info['calculator_parameters'] = data.get('calculator_parameters')
            del data['calculator_parameters']
        elif atoms.info.get('calculator_parameters'):
            pass
        else:
            raise ValueError('Calculation parameters missing.')

        for k, v in data.items():
            if isinstance(v, np.ndarray):
                fwio.array_to_list(v)
                data[k] = v

        encoding = fwio.atoms_to_encode(atoms)
        t0 = fireworks.PyTask(
            func='catkit.flow.fwio.encode_to_atoms',
            args=[encoding])

        t1 = fireworks.PyTask(
            func='catkit.flow.fwase.get_potential_energy',
            args=[calculator],
            stored_data_varname='trajectory')
        tasks = [t0, t1]

        if spec is None:
            spec = {'keys': keys, 'data': data}
        else:
            spec.update({'keys': keys, 'data': data})

        firework = fireworks.Firework(tasks, spec=spec)
        workflow = fireworks.Workflow([firework], name=workflow_name)
        self.launchpad.add_wf(workflow)


    def submit_relaxation_db(self, database, spec=None):
        """Submit each entry of an ASE database for relaxation.
        This requires that the calculation parameters be stored in
        the data under `data.calculator_parameters`.

        Parameters
        ----------
        database : str
            Path to ASE database to be looped over for relaxation.
        spec : dict
            Additional specification to be passed to Fireworks.
        """
        filename = database.split('/')[-1]

        db = connect(database)
        for d in db.select():
            self.submit_relaxation(
                d, calculation_name=filename, spec=spec)

    def bulk_relaxation(
            self,
            atoms,
            parameters,
            spec=None):
        """Run a relaxation of a given DB entry or atoms object. If a
        database object is used, the calculation will automatically store
        the keys and data for later retrieval.

        The entries uuid will also be stored and `data.calculator_parameters`
        will be used as the calculation parameters.

        Parameters
        ----------
        images : Atoms object
            Initial atoms to perform workflow on.
        parameters : dict
            Calculation parameters to use.
        workflow_name : str
            Name of the fireworks calculation to be used.
        spec : dict
            Additional fireworks specifications to pass to the database.
        """
        atoms.info['calculator_parameters'] = parameters

        encoding = fwio.atoms_to_encode(atoms)
        t0 = fireworks.PyTask(
            func='catkit.flow.fwio.encode_to_atoms',
            args=[encoding])

        t1 = fireworks.PyTask(
            func='catkit.flow.fwase.catflow_relaxation',
            stored_data_varname='trajectory')

        tasks = [t0, t1]

        if spec is None:
            spec = {}

        firework = fireworks.Firework(tasks, spec=spec)
        workflow = fireworks.Workflow([firework], name='bulk_relaxation')
        workflow_id = self.launchpad.add_wf(workflow)[-1]

        return workflow_id
