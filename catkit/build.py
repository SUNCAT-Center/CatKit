import catkit
import numpy as np
import ase.build
import ase


def bulk(name, crystalstructure=None, primitive=False, **kwargs):
    """Return the standard conventional cell of a bulk structure
    created using ASE. Accepts all keyword arguments for the ase
    bulk generator.

    Parameters
    ----------
    name : Atoms object | str
        Chemical symbol or symbols as in 'MgO' or 'NaCl'.
    crystalstructure : str
        Must be one of sc, fcc, bcc, hcp, diamond, zincblende,
        rocksalt, cesiumchloride, fluorite or wurtzite.
    primitive : bool
        Return the primitive unit cell instead of the conventional
        standard cell.

    Returns
    -------
    standardized_bulk : Gratoms object
        The conventional standard or primitive bulk structure.
    """
    if isinstance(name, str):
        atoms = ase.build.bulk(name, crystalstructure, **kwargs)
    else:
        atoms = name
    standardized_bulk = catkit.gen.symmetry.get_standardized_cell(
        atoms, primitive=primitive)

    return standardized_bulk


def surface(
        elements,
        size,
        crystal='fcc',
        miller=(1, 1, 1),
        termination=0,
        fixed=0,
        vacuum=10,
        orthogonal=False,
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
    orthogonal : bool
        Force the slab generator to produce the most orthogonal slab.

    Returns
    -------
    slab : Gratoms object
        Return a slab generated from the specified bulk structure.
    """
    if isinstance(elements, ase.Atoms):
        atoms = elements
    else:
        atoms = ase.build.bulk(elements, crystal, cubic=True, **kwargs)

    generator = catkit.gen.surface.SlabGenerator(
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
    elif len(size) == 3 and not orthogonal:
        size = size[:2]

    if orthogonal:
        catkit.gen.defaults['orthogonal'] = True
        if isinstance(size, (list, tuple)):
            size = np.prod(size[:2])

    slab = generator.get_slab(size=size, iterm=termination)

    return slab


def molecule(species, bond_index=None, vacuum=0):
    """Return list of enumerated gas-phase molecule structures based
    on species and topology.

    Parameters
    ----------
    species : str
        The chemical symbols to construct a molecule from.
    bond_index : int
        Construct the molecule as though it were adsorbed to a surface
        parallel to the z-axis. Will bond by the atom index given.
    vacuum : float
        Angstroms of vacuum to pad the molecules with.

    Returns
    -------
    images : list of Gratoms objects
        3D structures of the requested chemical species and topologies.
    """
    molecule_graphs = catkit.gen.molecules.get_topologies(species)

    images = []
    for atoms in molecule_graphs:
        atoms = catkit.gen.molecules.get_3D_positions(atoms, bond_index)
        images += [atoms]

    return images
