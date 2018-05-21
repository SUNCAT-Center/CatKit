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
        index=(1, 1, 1),
        fixed=0,
        vacuum=10,
        root=None,
        primitive=True,
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
    index : list (3,) or (4,)
        The miller index to cleave the surface structure from.
    fixed : int
        Number of layers to constrain.
    vacuum : float
        Angstroms of vacuum to add to the unit cell.
    root : int
        If not None, attempt to produce a root unit cell with
        a primitive lattice length multiple this root.
    primitive : bool
        Perform an spglib reduction of the slabs unit cell.

    Returns
    -------
    slab : object
        Return a slab generated from the specified bulk structure.
    """
    if isinstance(elements, Atoms):
        atoms = elements
    else:
        atoms = bulk(elements, crystal, cubic=True, **kwargs)

    gen = SlabGenerator(
        atoms,
        miller_index=index,
        layers=size[2],
        fixed=fixed,
        vacuum=vacuum)

    slab = gen.get_slab(size=size, root=root, primitive=primitive)
    surface_atoms = gen.get_voronoi_surface_atoms(slab, attach_graph=True)
    slab.set_surface_atoms(surface_atoms[0])

    return slab


def molecule(
        species,
        topology=None,
        vacuum=0):
    """Return gas-phase molecule structures based on species and
    topology.

    Parameters
    ----------
    species : str
        The chemical symbols to construct a molecule from.
    topology : int, str, or slice
        The indices for the distinct topology produced by the generator.
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
        molecule_graphs = get_topologies(species)[_slice]

    images = []
    for atoms in molecule_graphs:
        branches = nx.bfs_successors(atoms.graph, 0)

        root = None
        for i, branch in enumerate(branches):
            _branch_molecule(atoms, branch, base_root=root)
            root = 0

        if vacuum:
            atoms.center(vacuum=vacuum)
        images += [atoms]

    return images
