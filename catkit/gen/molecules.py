import catkit
import networkx as nx
import itertools
import numpy as np


def bin_hydrogen(hydrogens=1, bins=1):
    """Recursive function for determining distributions of
    hydrogens across bins.
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
    atoms += catkit.Gratoms('H{}'.format(sum(bins)))
    atoms.graph.add_edges_from(edges)

    return atoms


def get_topologies(symbols, saturate=False):
    """Return the possible topologies of a given chemical species.

    Parameters
    ----------
    symbols : str
        Atomic symbols to construct the topologies from.
    saturate : bool
        Saturate the molecule with hydrogen based on the
        default.radicals set.

    Returns
    -------
    molecules : list (N,)
        Gratoms objects with unique connectivity matrix attached.
        No 3D positions will be provided for these structures.
    """
    num, cnt = catkit.gen.utils.get_atomic_numbers(symbols, True)
    mcnt = cnt[num != 1]
    mnum = num[num != 1]

    if cnt[num == 1]:
        hcnt = cnt[num == 1][0]
    else:
        hcnt = 0

    elements = np.repeat(mnum, mcnt)
    max_degree = catkit.gen.defaults.get('radicals')[elements]
    n = mcnt.sum()

    hmax = int(max_degree.sum() - (n - 1) * 2)
    if hcnt > hmax:
        hcnt = hmax

    if saturate:
        hcnt = hmax

    if n == 1:
        atoms = catkit.Gratoms(elements, cell=[1, 1, 1])
        hatoms = hydrogenate(atoms, np.array([hcnt]))
        return [hatoms]
    elif n == 0:
        hatoms = catkit.Gratoms('H{}'.format(hcnt))
        if hcnt == 2:
            hatoms.graph.add_edge(0, 1, bonds=1)
        return [hatoms]

    ln = np.arange(n).sum()
    il = np.tril_indices(n, -1)

    backbones, molecules = [], []
    combos = itertools.combinations(np.arange(ln), n - 1)
    for c in combos:
        # Construct the connectivity matrix
        ltm = np.zeros(ln)
        ltm[np.atleast_2d(c)] = 1

        connectivity = np.zeros((n, n))
        connectivity[il] = ltm
        connectivity = np.maximum(connectivity, connectivity.T)

        degree = connectivity.sum(axis=0)

        # Not fully connected (subgraph)
        if np.any(degree == 0) or not \
           nx.is_connected(nx.from_numpy_array(connectivity)):
            continue

        # Overbonded atoms.
        remaining_bonds = (max_degree - degree).astype(int)
        if np.any(remaining_bonds < 0):
            continue

        atoms = catkit.Gratoms(
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


def get_3D_positions(atoms, bond_index=None):
    """Return an estimation of the 3D structure of a Gratoms object
    based on its graph.

    WARNING: This function operates on the atoms object in-place.

    Parameters
    ----------
    atoms : Gratoms object
        Structure with connectivity matrix to provide a 3D structure.
    bond_index : int
        Index of the atoms to consider as the origin of a surface
        bonding site.

    Returns
    -------
    atoms : Gratoms object
        Structure with updated 3D positions.
    """
    branches = nx.dfs_successors(atoms.graph, bond_index)

    complete = []
    for i, branch in enumerate(branches.items()):
        root, nodes = branch

        if len(nodes) == 0:
            continue

        c0 = atoms[root].position
        if i == 0:
            basis = catkit.gen.utils.get_basis_vectors([c0, [0, 0, -1]])
        else:
            bond_index = None
            for j, base_root in enumerate(complete):
                if root in branches[base_root]:
                    c1 = atoms[base_root].position
                    # Flip the basis for every alternate step down the chain.
                    basis = catkit.gen.utils.get_basis_vectors([c0, c1])
                    if (i - j) % 2 != 0:
                        basis[2] *= -1
                    break
        complete.insert(0, root)

        positions = _branch_molecule(atoms, branch, basis, bond_index)
        atoms.positions[nodes] = positions

    return atoms


def _branch_molecule(
        atoms,
        branch,
        basis=None,
        adsorption=None):
    """Return the positions of a Gratoms object for a segment of its
    attached graph. This function is mean to be iterated over by a depth
    first search form NetworkX.

    Parameters
    ----------
    atoms : Gratoms object
        Gratoms object to be iterated over. Will have its positions
        altered in-place.
    branch : tuple (1, [N,])
        A single entry from the output of nx.bfs_successors. The
        first entry is the root node and the following list is the
        nodes branching from the root.
    basis : ndarray (3, 3)
        The basis vectors to use for this step of the branching.
    adsorption : bool
        If True, will construct the molecule as though there is a
        surface normal to the negative z-axis. Must be None for
        all but the first index in the depth first search.

    Returns
    -------
    positions : ndarray (N, 3)
        Estimated positions for the branch nodes.
    """
    root, nodes = branch
    root_position = atoms[root].position

    radii = catkit.gen.defaults.get('radii')
    atomic_numbers = atoms.numbers[[root] + nodes]
    atomic_radii = radii[atomic_numbers]
    dist = (atomic_radii[0] + atomic_radii[1:])[:, None]

    angles = np.array([109.47, 109.47, 109.47, 0])
    dihedral = np.array([0, 120, -120, 0])

    if adsorption:
        # Move adsorption structures away from surface
        angles += 15

    # Tetrahedral bond arrangement by default
    n = len(nodes)
    if n == 1:
        # Linear bond arrangement
        angles[0] = 180
    elif n == 2:
        # Trigonal-planer bond arrangement
        angles[:2] = 120
        dihedral[1] = 180

    # Place the atoms of this segment of the branch
    if basis is None:
        basis = catkit.gen.utils.get_basis_vectors(
            [root_position, [0, 0, -1]])
    basis = np.repeat(basis[None, :, :], len(dist), axis=0)

    ang = np.deg2rad(angles)[:len(dist), None]
    tor = np.deg2rad(dihedral)[:len(dist), None]

    basis[:, 1] *= -np.sin(tor)
    basis[:, 2] *= np.cos(tor)

    vectors = basis[:, 1] + basis[:, 2]
    vectors *= dist * np.sin(ang)
    basis[:, 0] *= dist * np.cos(ang)

    positions = vectors + root_position - basis[:, 0]

    return positions
