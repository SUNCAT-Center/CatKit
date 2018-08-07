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
            password=None,
            calculator='espresso.Espresso'
    ):
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
        self.calculator = calculator

    def submit_relaxation(
            self,
            image,
            calculation_name,
            parameters=None,
            spec=None):
        """Run a relaxation of a given DB entry or atoms object. If a
        database object is used, the calculation will automatically store
        the keys and data for later retrieval.

        The entries uuid will also be stored and `data.calculator_parameters`
        will be used as the calculation parameters.

        Parameter:
        ----------
        images : Atoms object | AtomsRow object
            ASE database entry or atoms object to relax.
        calculation_name : str
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

        if parameters is not None:
            atoms.info = parameters
        elif data.get('calculator_parameters'):
            atoms.info = data.get('calculator_parameters')
        else:
            raise ValueError('Calculation parameters missing.')

        # These will be stored in the atoms info.
        del data['calculator_parameters']

        for k, v in data.items():
            if isinstance(v, np.ndarray):
                fwio.array_to_list(v)
                data[k] = v

        encoding = fwio.atoms_to_encode(atoms)
        t0 = fireworks.PyTask(
            func='catkit.flow.fwio.encode_to_atoms',
            args=[encoding])

        t1 = fireworks.PyTask(
            func='catkit.flow.fw_ase.get_potential_energy',
            args=[self.calculator],
            stored_data_varname='trajectory')

        if spec is None:
            spec = {'keys': keys, 'data': data}
        else:
            spec.update({'keys': keys, 'data': data})

        firework = fireworks.Firework(
            [t0, t1], spec=spec,
            name=calculation_name)

        workflow = fireworks.Workflow([firework])
        self.launchpad.add_wf(workflow)

    def submit_relaxation_db(self, database, spec=None):
        """Submit each entry of an ASE database for relaxation.
        This requires that the calculation parameters be stored in
        the data under `data.calculator_parameters`.

        Parameter:
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
