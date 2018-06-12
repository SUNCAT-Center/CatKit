from .gen.surface import SlabGenerator
from .gen.molecules import get_topologies
from .gen.geometry import _branch_molecule
from .gen import utils
from ase.build import bulk
import networkx as nx
from ase import Atoms


def surface(
        elements,
        size,
        crystal='fcc',
        miller=(1, 1, 1),
        termination=0,
        fixed=0,
        vacuum=10,
        **kwargs):
    """A helper function to return the surface associated with a
    given set of input parameters to the general surface generator.

    Parameters
    ----------
    elements : str or object
        The atomic symbol to be passed to the as bulk builder function
        or an atoms object representing the bulk structure to use.
    size : list (3,)
        Number of time to expand the x, y, and z primitive cell.
    crystal : str
        The bulk crystal structure to pass to the ase bulk builder.
    miller : list (3,) or (4,)
        The miller index to cleave the surface structure from. If 4 values
        are used, assume Miller-Bravis convention.
    termination : int
        The index associated with a specific slab termination.
    fixed : int
        Number of layers to constrain.
    vacuum : float
        Angstroms of vacuum to add to the unit cell.

    Returns
    -------
    slab : Gratoms object
        Return a slab generated from the specified bulk structure.
    """
    if isinstance(elements, Atoms):
        atoms = elements
    else:
        atoms = bulk(elements, crystal, cubic=True, **kwargs)

    gen = SlabGenerator(
        bulk=atoms,
        miller_index=miller,
        layers=size[-1],
        vacuum=vacuum,
        fixed=fixed,
        layer_type=kwargs.get('layer_type', 'trim'),
        attach_graph=kwargs.get('attach_graph', True),
        standardize_bulk=kwargs.get('standardize_bulk', True),
        tol=kwargs.get('tol', 1e-8)
    )

    if len(size) == 2:
        size = size[0]
    elif len(size) == 3:
        size = size[:2]

    slab = gen.get_slab(size=size, iterm=termination)

    return slab


def add_adsorbate(atoms):
    """Add an adsorbate to a surface."""

    return atoms


def molecule(
        species,
        topology=None,
        adsorption=False,
        vacuum=0):
    """Return gas-phase molecule structures based on species and
    topology.

    Parameters
    ----------
    species : str
        The chemical symbols to construct a molecule from.
    topology : int, str, or slice
        The indices for the distinct topology produced by the generator.
    adsorption : bool
        Construct the molecule as though it were adsorbed to a surface
        parallel to the z-axis.
    vacuum : float
        Angstroms of vacuum to pad the molecule with.

    Returns
    -------
    images : list of objects
        3D structures of the requested chemical species and topologies.
    """
    molecule_graphs = get_topologies(species)

    if len(molecule_graphs) > 1:
        _slice = utils.parse_slice(topology)
        molecule_graphs = molecule_graphs[_slice]

    images = []
    for atoms in molecule_graphs:
        branches = nx.bfs_successors(atoms.graph, 0)

        root = None
        for i, branch in enumerate(branches):
            _branch_molecule(atoms, branch, root, adsorption)
            root = 0

        if vacuum:
            atoms.center(vacuum=vacuum)
        images += [atoms]

    return images
