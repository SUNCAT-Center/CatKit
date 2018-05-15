from catgen.surface import SlabGenerator
from ase.build import bulk
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

    Parameters:
    -----------
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

    Returns:
    --------
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

    return slab
