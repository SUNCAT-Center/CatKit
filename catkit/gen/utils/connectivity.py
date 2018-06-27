from .. import defaults
from . import coordinates
from scipy.spatial import Voronoi
import warnings
import numpy as np


def get_voronoi_neighbors(atoms, r=4):
    """Return the connectivity matrix from the Voronoi
    method. Multi-bonding occurs through periodic boundary conditions.

    Parameters
    ----------
    atoms : atoms object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.
    r : float
        Radius of the spheres to expand around each atom.

    Returns
    -------
    connectivity : ndarray (n, n)
        Number of edges formed between atoms in a system.
    """
    index, coords, offsets = coordinates.expand_cell(atoms, r)

    L = int(len(offsets) / 2)
    origional_indices = np.arange(L * len(atoms), (L + 1) * len(atoms))

    voronoi = Voronoi(coords, qhull_options='QbB Qc Qs')
    points = voronoi.ridge_points

    connectivity = np.zeros((len(atoms), len(atoms)))
    for i, n in enumerate(origional_indices):
        p = points[np.where(points == n)[0]]
        edges = np.sort(index[p])

        if not edges.size:
            warnings.warn(
                ("scipy.spatial.Voronoi returned an atom which has "
                 "no neighbors. This may result in incorrect connectivity."))
            continue

        unique_edge, edge_counts = np.unique(
            edges,
            return_counts=True,
            axis=0)

        for j, edge in enumerate(unique_edge):
            u, v = edge
            connectivity[u][v] += edge_counts[j]
            connectivity[v][u] += edge_counts[j]

    connectivity /= 2

    return connectivity.astype(int)


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
    cov_radii = defaults.get('covalent_radii')
    numbers = atoms.numbers
    index, coords = coordinates.expand_cell(atoms, 4)[:2]

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
