from catkit import Gratoms
from .utils import get_atomic_numbers
from itertools import combinations
import numpy as np

radicals = np.ones(92)
radicals[[6, 7, 8, 9, 15, 16]] = [4, 3, 2, 1, 3, 2]


def bin_hydrogen(hydrogens=1, bins=1):
    """Recursive function for determining distributions of
    hydrogens across bins
    """
    if bins == 1:
        yield [hydrogens]

    elif hydrogens == 0:
        yield [0] * bins

    else:
        for i in range(hydrogens + 1):
            for j in bin_hydrogen(hydrogens - i, 1):
                for k in bin_hydrogen(i, bins - 1):
                    yield j + k


def hydrogenate(atoms, bins, copy=True):
    """Add hydrogens to a gratoms object via provided bins"""
    h_index = len(atoms)

    edges = []
    for i, j in enumerate(bins):
        for _ in range(j):
            edges += [(i, h_index)]
            h_index += 1

    if copy:
        atoms = atoms.copy()
    atoms += Gratoms('H{}'.format(sum(bins)))
    atoms.graph.add_edges_from(edges)

    return atoms


def get_topologies(chemistry, saturate=False):
    """Return the possible topologies of a given chemical species
    """
    num, cnt = get_atomic_numbers(chemistry, True)
    mcnt = cnt[num != 1]
    mnum = num[num != 1]

    if cnt[num == 1]:
        hcnt = cnt[num == 1][0]
    else:
        hcnt = 0

    elements = np.repeat(mnum, mcnt)
    max_degree = radicals[elements]
    n = mcnt.sum()

    hmax = int(max_degree.sum() - (n - 1) * 2)
    if hcnt > hmax:
        hcnt = hmax

    if saturate:
        hcnt = hmax

    if n == 1:
        atoms = Gratoms(elements, cell=[1, 1, 1])
        hatoms = hydrogenate(atoms, np.array([hcnt]))
        return [hatoms]
    elif n == 0:
        hatoms = Gratoms('H{}'.format(hcnt))
        if hcnt == 2:
            hatoms.graph.add_edge(0, 1, bonds=1)
        return [hatoms]

    ln = np.arange(n).sum()
    il = np.tril_indices(n, -1)

    backbones, molecules = [], []
    for c in combinations(np.arange(ln), n - 1):
        # Construct the connectivity matrix
        ltm = np.zeros(ln)
        ltm[[c]] = 1

        connectivity = np.zeros((n, n))
        connectivity[il] = ltm
        connectivity = np.maximum(connectivity, connectivity.T)

        degree = connectivity.sum(axis=0)

        # Not fully connected (cyclical subgraph)
        if np.any(degree == 0):
            continue

        # Overbonded atoms.
        remaining_bonds = (max_degree - degree).astype(int)
        if np.any(remaining_bonds < 0):
            continue

        atoms = Gratoms(
            numbers=elements,
            edges=connectivity,
            cell=[1, 1, 1])

        isomorph = False
        for G0 in backbones:
            if atoms.is_isomorph(G0):
                isomorph = True
                break

        if not isomorph:
            backbones += [atoms]

            # The backbone is saturated, do not enumerate
            if hcnt == hmax:
                hatoms = hydrogenate(atoms, remaining_bonds)
                molecules += [hatoms]
                continue

            # Enumerate hydrogens across backbone
            for bins in bin_hydrogen(hcnt, n):
                if not np.all(bins <= remaining_bonds):
                    continue

                hatoms = hydrogenate(atoms, bins)

                isomorph = False
                for G0 in molecules:
                    if hatoms.is_isomorph(G0):
                        isomorph = True
                        break

                if not isomorph:
                    molecules += [hatoms]

    return molecules
