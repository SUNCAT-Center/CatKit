from . import operations
from .. import Gratoms
from ..gen import utils
import pkg_resources
import json
import numpy as np
import ase


class Fingerprinter():
    """Parent class for all fingerprint generators."""

    def __init__(self, images=None):
        """Store the properties of the databases being loaded as
        as needed.

        Parameters
        ----------
        images : list
            Atoms objects to generate fingerprints for.
        """
        if isinstance(images, ase.Atoms):
            images = [images]
        self._images = images

        self._data = None
        self._slab_d_data = None
        self._bulk_d_data = None
        self._mendeleev_data = None

    def _get_data(self):
        """A lazy evaluator for loading general data."""
        data = self._data

        if data is None:
            path = pkg_resources.resource_filename(
                'catkit', 'data/properties.json')

            with open(path, 'r') as f:
                data = json.load(f)

            self._data = data

        return data

    def _get_atoms_parameters(self, atoms, parameters):
        """Return the parameters associated with a specific atoms
        object and seed parameters set.

        Parameters
        ----------
        atoms : Atoms object
            Structure to return parameters for, contains n atoms.
        parameters : list of str (M,)
            Seed parameters to use for atom specific parameter
            generation.

        Returns
        -------
        atoms_parameters : ndarray (M, N)
            General chemical properties specific to this atoms object
            and list of seed parameters.
        """
        data = self._get_data()

        num = atoms.numbers
        atoms_parameters = np.zeros((len(parameters), len(num)))
        for i, k in enumerate(parameters):
            atoms_parameters[i] = np.asarray(data[k])[num]

        return atoms_parameters

    def get_fp(self, parameters, operation_list):
        """Return the fingerprints for a list of images of single atoms
        object for the given parameters. Convolutions will be performed for
        specific operations if the parameters is provided as a list of lists.

        Parameters
        ----------
        parameters : list of str | list of lists (M,)
            Names of seeding parameters available in the parameters
            database. If a list lists is provided, the number of lists
            must be equal to the number of operations.
        operation_list : list of func | str (M,)
            A list of operation functions to produce the fingerprints from.
            The names (str) of operations functions can also be used.

        Returns
        -------
        fingerprints : ndarray (N, X)
            Fingerprints for the images produced from the provided
            seed parameters.
        """
        fingerprints = [[] for _ in self._images]
        for i, atoms in enumerate(self._images):

            method = None
            if np.all(atoms.pbc):
                method = 'voronoi'

            if isinstance(atoms, Gratoms):
                connectivity = atoms.connectivity
            else:
                connectivity = get_connectivity(atoms, method=method)

            for j, operation in enumerate(operation_list):
                kwargs = {}
                if isinstance(operation, list):
                    operation, kwargs = operation

                if isinstance(operation, str):
                    operation = getattr(operations, operation)

                if all(isinstance(pl, list) for pl in parameters):
                    atoms_parameters = self._get_atoms_parameters(
                        atoms, parameters[j])
                else:
                    atoms_parameters = self._get_atoms_parameters(
                        atoms, parameters)

                fingerprint = _generate_fingerprint(
                    operation,
                    atoms,
                    atoms_parameters,
                    connectivity,
                    **kwargs)
                fingerprints[i] += [fingerprint]
        fingerprints = np.block(fingerprints)

        return fingerprints


def _generate_fingerprint(
        operation,
        atoms,
        atoms_parameters,
        connectivity,
        **kwargs):
    """Return the parameters of a convolutions of atomic properties
    based on connectivity in a periodic unit cell.

    Parameters
    ----------
    operations : list of functions
        A list of operation functions to produce the fingerprints from.
    atoms : object
        Atomic structure to generate fingerprint for.
    atoms_parameters : ndarray (N, M)
        General chemical properties specific to this atoms object
        and list of seed parameters.
    connectivity : ndarray (N, N)
        Estimated connectivity matrix where n is the number
        of atoms in the atoms-object.

    Returns
    -------
    fingerprint : ndarray (M,)
        Fingerprint produced from a given operation and seed parameters.
    """
    fingerprint = operation(
        atoms=atoms,
        atoms_parameters=atoms_parameters,
        connectivity=connectivity,
        **kwargs)

    return fingerprint


def get_connectivity(atoms, method=None):
    """Returns an estimate of the connectivity matrix
    for a given atoms-object from CatGen.

    Parameters
    ----------
    atoms : object
        Molecular structure with out without adsorbates.
    method : str (None or 'voronoi')
        Method for estimating the connectivity matrix:

        None - standard cutoff radius method.
        voronoi - best suited for bulk characterization.

    Returns
    -------
    connectivity : ndarray (N, N)
        Estimated connectivity matrix where n is the number
        of atoms in the atoms-object.
    """
    if method == 'voronoi':
        connectivity = utils.get_voronoi_neighbors(atoms)
    else:
        connectivity = utils.get_cutoff_neighbors(atoms)

    return connectivity
