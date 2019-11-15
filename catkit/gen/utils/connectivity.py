import catkit
from . import coordinates
import numpy as np
import scipy
import warnings


def get_voronoi_neighbors(atoms, cutoff=5.0, return_distances=False):
    """Return the connectivity matrix from the Voronoi
    method. Multi-bonding occurs through periodic boundary conditions.

    Parameters
    ----------
    atoms : atoms object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.
    cutoff : float
        Radius of maximum atomic bond distance to consider.

    Returns
    -------
    connectivity : ndarray (n, n)
        Number of edges formed between atoms in a system.
    """
    index, coords, offsets = coordinates.expand_cell(atoms, cutoff=cutoff)

    L = int(len(offsets) / 2)
    origional_indices = np.arange(L * len(atoms), (L + 1) * len(atoms))

    voronoi = scipy.spatial.Voronoi(coords, qhull_options='QbB Qc Qs')
    points = voronoi.ridge_points

    connectivity = np.zeros((len(atoms), len(atoms)))
    distances = []
    distance_indices = []
    for i, n in enumerate(origional_indices):
        p = points[np.where(points == n)[0]]
        dist = np.linalg.norm(np.diff(coords[p], axis=1), axis=-1)[:, 0]
        edges = np.sort(index[p])

        if not edges.size:
            warnings.warn(
                ("scipy.spatial.Voronoi returned an atom which has "
                 "no neighbors. This may result in incorrect connectivity."))
            continue

        unique_edge = np.unique(edges, axis=0)

        for j, edge in enumerate(unique_edge):
            indices = np.where(np.all(edge == edges, axis=1))[0]

            d = dist[indices][np.where(dist[indices] < cutoff)[0]]
            count = len(d)
            if count == 0:
                continue

            u, v = edge

            distance_indices += [sorted([u, v])]
            distances += [sorted(d)]

            connectivity[u][v] += count
            connectivity[v][u] += count

    connectivity /= 2

    if not return_distances:
        return connectivity.astype(int)

    distance_indices, unique_idx_idx = \
        np.unique(distance_indices, axis=0, return_index=True)

    distances = [distances[i] for i in unique_idx_idx]

    pair_distances = {'indices': distance_indices.tolist(),
                      'distances': distances}

    if return_distances:
        return connectivity.astype(int), pair_distances


def get_cutoff_neighbors(atoms, cutoff=None, atol=1e-8):
    """Return the connectivity matrix from a simple radial cutoff.
    Multi-bonding occurs through periodic boundary conditions.

    Parameters
    ----------
    atoms : atoms object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.
    cutoff : float
        Cutoff radius to use when determining neighbors.
    atol: float
        Absolute tolerance to use when computing distances.

    Returns
    -------
    connectivity : ndarray (n, n)
        Number of edges formed between atoms in a system.
    """
    cov_radii = catkit.gen.defaults.get('radii')
    numbers = atoms.numbers
    index, coords = coordinates.expand_cell(atoms)[:2]

    if cutoff is None:
        radii = cov_radii[numbers]
        radii = np.repeat(radii[:, None], len(index) / len(radii), axis=1)
        radii = radii.T.flatten()
    else:
        radii = np.ones_like(index) * cutoff

    connectivity = np.zeros((len(atoms), len(atoms)), dtype=int)
    for i, center in enumerate(atoms.positions):
        norm = np.linalg.norm(coords - center, axis=1)
        boundry_distance = norm - radii
        if cutoff is None:
            boundry_distance -= cov_radii[numbers[i]]
        neighbors = index[np.where(boundry_distance < atol)]

        unique, counts = np.unique(neighbors, return_counts=True)
        connectivity[i, unique] = counts
        # Remove self-interaction
        connectivity[i, i] -= 1

    return connectivity
