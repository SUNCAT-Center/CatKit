from .gen.surface import SlabGenerator
from .gen.surface import get_unique_indices
from .gen import utils
import numpy as np


def surfaces(
        bulk,
        width,
        miller_indices=(1, 1, 1),
        terminations=None,
        sizes=None,
        vacuum=10,
        fixed=0,
        layer_type='angs',
        **kwargs):
    """Return a list of enumerated surfaces based on symmetry properties of
    interest to the user. Any bulk structure provided will be standardized.

    This function will take additional keyword arguments for the
    :meth:`catkit.gen.surface.SlabGenerator` Class.

    Parameters
    ----------
    bulk : str | Atoms
        The atomic symbol to be passed to the as bulk builder function
        or an atoms object representing the bulk structure to use.
    width : float
        Minimum width of the slab in angstroms before trimming. Imposing
        symmetry requirements will reduce the width.
    miller_indices : int | list (3,) | list of list (n, 3)
        List of the miller indices to enumerate slabs for. If an integer is
        provided, the value is treated as the maximum miller index to consider
        for an enumeration of all possible unique miller indices.
    terminations : int | array_like
        Return the terminations associated with the provided indices. If -1,
        all possible terminations are enumerated.
    sizes : None | int | array_like (n,)
        Enumerate all surface sizes in the provided list. Sizes are integers
        which represent multiples of the smallest possible surface area.
        If None, return slabs with the smallest possible surface area. If an
        integer, enumerate all sizes up to that multiple.
    vacuum : float
        Angstroms of vacuum to add to the unit cell.
    fixed : int
        Number of layers to constrain.
    layer_type : 'angs' | 'trim' | 'stoich' | 'sym'
        Method of slab layering to perform. See also:
        :meth:`catkit.gen.surface.SlabGenerator`

    Returns
    -------
    slabs : list of Gratoms objects
        Return a list of enumerated slab structures.
    """
    standardized_bulk = utils.get_spglib_cell(bulk, tol=5e-3)

    if isinstance(miller_indices, int):
        miller_indices = get_unique_indices(standardized_bulk, miller_indices)
    elif isinstance(miller_indices, (list, np.ndarray)):
        miller_indices = np.atleast_2d(miller_indices)

    if sizes is None:
        sizes = np.ones(1)
    elif isinstance(sizes, int):
        sizes = np.arange(sizes)

    slabs = []
    for miller in miller_indices:
        gen = SlabGenerator(
            bulk=standardized_bulk,
            miller_index=miller,
            layers=width,
            vacuum=vacuum,
            fixed=fixed,
            layer_type=layer_type,
            **kwargs)

        if terminations is None:
            iterms = np.zeros(1)
        elif terminations == -1:
            zshifts = gen.get_unique_terminations()
            iterms = np.arange(len(zshifts))
        else:
            iterms = terminations

        for i in iterms:
            for size in sizes:
                slab = gen.get_slab(size=int(size), iterm=i)
                slab.info['miller'] = miller
                slabs += [slab]

    return slabs
