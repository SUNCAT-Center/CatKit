from . import utils
from ase import Atoms
from scipy.spatial import Delaunay
from scipy.linalg import circulant
from itertools import product
from numpy.linalg import pinv, norm
import numpy as np


def get_adsorption_sites(
        slab,
        surface_sites,
        symmetry_reduced=True,
        vectors=False,
        return_proxy=False,
        tol=1e-5,
):
    """ Get the adsorption sites of a slab as defined by surface
    atom symmetries.

    Parameters:
      slab: ASE atoms-object
        The slab to find adsorption sites for.
      surface_sites: ndarray (n,)
        surface sites of the provided slab.
      symmetry_reduced: int
        Whether to return the symmetrically unique sites only.
      vectors: bool
        Whether to compute the adsorption vectors.
      return_proxy: bool
        Whether to return the proxy_slab for adsorption work.

    Returns: dict of 3 lists
      Dictionary of top, bridge, hollow, and 4-fold sites containing
      positions, points, and neighbor lists. If adsorption vectors
      are requested, the third list is replaced.
    """

    sites, vslab = find_adsorption_sites(
        slab=slab,
        surface_sites=surface_sites,
        trim=0.5,
        tol=tol,
    )

    sites = get_reduced_sites(
        sites=sites,
        slab=slab,
        tol=tol,
    )

    if symmetry_reduced:
        sites = get_symmetric_sites(
            sites=sites,
            slab=slab,
            tol=tol,
        )

    if vectors:
        sites = _get_adsorption_vectors(vslab, sites)

    if return_proxy:
        return sites, vslab

    return sites


def find_adsorption_sites(
        slab,
        surface_sites,
        trim=0.5,
        tol=1e-5,
):
    """ Find all bridge and hollow sites (3-fold and 4-fold) given an
    input slab based Delaunay triangulation of surface atoms of a
    super-cell.

    Parameters:
      slab: ASE atoms-object
        The slab to find adsorption sites for.
      surface_sites: ndarray (n,)
        surface sites of the provided slab.
      trim: float
        Percentage of fractional coordinates to remove.
      tol: float
        Absolute tolerance for floating point errors.

    Returns: dict of 3 lists
      Dictionary of top, bridge, hollow, and 4-fold sites containing
      positions, points, and neighbors of the super-cell.
    """

    # Top sites projected into expanded unit cell
    index, coords, offsets = utils.expand_cell(slab, r=8)

    if isinstance(surface_sites, list):
        surface_sites = np.array(surface_sites)

    # Some magic to get the top sites indices
    n = len(slab)
    rtop = surface_sites[:, None].repeat(len(offsets), axis=1)
    rtop = (rtop + np.arange(0, len(coords), n)).T.flatten()

    top_coords = coords[rtop]
    numbers = slab.get_atomic_numbers().tolist() * len(offsets)

    vslab = Atoms(
        positions=top_coords,
        symbols=np.array(numbers)[rtop],
        cell=slab.cell,
        tags=rtop,
    )

    # Dict of 2 lists (positions, site_id)
    sites = {
        'top': [
            top_coords,
            np.arange(rtop.shape[0]),
            [[] for _ in range(rtop.shape[0])]
        ],
        'bridge': [[], [], []],
        'hollow': [[], [], []],
        '4fold': [[], [], []],
    }

    dt = Delaunay(sites['top'][0][:, :2])
    neighbors = dt.neighbors
    simplices = dt.simplices

    for i, corners in enumerate(simplices):

        cir = circulant(corners)
        edges = cir[:, 1:]

        # Inner angle of each triangle corner
        vec = sites['top'][0][edges.T] - sites['top'][0][corners]
        uvec = vec.T / norm(vec, axis=2).T
        angles = np.sum(uvec.T[0] * uvec.T[1], axis=1)

        # Angle types
        right = np.isclose(angles, 0)
        obtuse = (angles < -tol)

        rh_corner = corners[right]
        edge_neighbors = neighbors[i]

        if obtuse.any():
            # Assumption: All simplices with obtuse angles
            # are irrelevant boundaries.
            continue

        bridge = np.sum(sites['top'][0][edges], axis=1) / 2.0

        # Looping through corners allows for elimination of
        # redundant points, identification of 4-fold hollows,
        # and collection of bridge neighbors.
        for j, c in enumerate(corners):

            edge = sorted(edges[j])

            if edge in sites['bridge'][1]:
                continue

            # Get the bridge neighbors (for adsorption vector)
            neighbor_simplex = simplices[edge_neighbors[j]]
            oc = list(set(neighbor_simplex) - set(edge))[0]

            # Right angles potentially indicate 4-fold hollow
            potential_hollow = edge + sorted([c, oc])
            # print(potential_hollow)
            if c in rh_corner:

                if potential_hollow in sites['4fold'][1]:
                    continue

                # Assumption: If not 4-fold, this suggests
                # no hollow OR bridge site is present.

                ovec = sites['top'][0][edge] - sites['top'][0][oc]
                ouvec = ovec / norm(ovec)
                oangle = np.dot(*ouvec)
                oright = np.isclose(oangle, 0)
                if oright:
                    sites['4fold'][0] += [bridge[j]]
                    sites['4fold'][1] += [potential_hollow]
                    sites['top'][2][c] += [oc]
            else:
                sites['bridge'][0] += [bridge[j]]
                sites['bridge'][1] += [edge]
                sites['bridge'][2] += [[c, oc]]

            sites['top'][2][edge[0]] += [edge[1]]
            sites['top'][2][edge[1]] += [edge[0]]

        if not right.any():
            hollow = np.average(sites['top'][0][corners], axis=0)
            sites['hollow'][0] += [hollow]
            sites['hollow'][1] += [corners.tolist()]

    # For collecting missed bridge neighbors
    for s in sites['4fold'][1]:

        for edge in product(s[:2], s[2:]):
            edge = sorted(edge)
            i = sites['bridge'][1].index(edge)
            n, m = sites['bridge'][1][i], sites['bridge'][2][i]
            nn = (set(s) - set(n + m))
            sites['bridge'][2][i] += [list(nn)[0]]

    xlim = offsets.T[0].max()
    ylim = offsets.T[1].max()

    # Convert lists to arrays
    for k, v in sites.items():
        positions, points, _ = v

        if len(positions) == 0:
            continue

        frac_coords = np.dot(positions, pinv(slab.cell))

        screen = (
            frac_coords[:, 0] < xlim + 1 - (trim * xlim)) & \
            (frac_coords[:, 0] > -xlim + (trim * xlim)) & \
            (frac_coords[:, 1] < ylim + 1 - (trim * ylim)) & \
            (frac_coords[:, 1] > -ylim + (trim * ylim))

        sites[k][0] = np.array(positions)[screen]
        sites[k][1] = np.array(points)[screen]

        if k in ['top', 'bridge']:
            sites[k][2] = np.array(sites[k][2])[screen]

    return sites, vslab


def get_reduced_sites(
        sites,
        slab,
        tol=1e-5,
):
    """ Reduce overlapping points via fractional coordinates. Intended
    for use after finding super-cell sites.

    Parameters:
      sites: dict
        Dictionary of cartesian coordinates to provide reduced sites from.
        Must have the form {'site': [[positions], [points], [neighbors]]}.
      slab: ASE atoms-object
        slab to determine symmetry operations from.
      tol: float
        Absolute tolerance for floating point errors.

    Returns: dict of 3 lists
      Dictionary sites containing positions, points, and neighbor lists.
    """

    for k, v in sites.items():
        positions, points, nnneighbors = v

        if len(positions) == 0:
            continue

        frac_coords = np.dot(positions, pinv(slab.cell))

        non_unique, unique = [], []
        for i, xyz in enumerate(frac_coords):
            if i in non_unique:
                continue

            non_unique += matching_sites(xyz, frac_coords).tolist()
            unique += [i]

        unique = np.array(unique)
        unique_positions = frac_coords[unique]

        shift = [tol, tol, 0]
        unique_positions += shift
        unique_positions %= 1
        unique_positions -= shift

        srt = np.lexsort((
            unique_positions[:, 0],
            unique_positions[:, 1]
        ))

        sites[k][0] = np.dot(unique_positions, slab.cell)[srt]
        sites[k][1] = points[unique][srt]

        if isinstance(nnneighbors, np.ndarray):
            sites[k][2] = nnneighbors[unique][srt]

    return sites


def get_symmetric_sites(
        sites,
        slab,
        tol=1e-5
):
    """ Determine the symmetrically unique adsorption sites
    from a dictionary of possible sites for a given slab.

    Parameters:
      sites: dict
        Dictionary of Cartesian coordinates to provide reduced sites from.
        Must have the form {'site': [[positions], [points], [neighbors]]}.
      slab: ASE atoms-object
        slab to determine symmetry operations from.
      tol: float
        Absolute tolerance for floating point errors.

    Returns: dict of 3 lists
      Dictionary sites containing positions, points, and neighbor lists.
    """

    symmetry = utils.get_symmetry(slab, tol=tol)
    rotations = symmetry['rotations']
    translations = symmetry['translations']

    for k, v in sites.items():
        positions, points, nnneighbors = v

        if len(positions) == 0:
            continue

        frac_coords = np.dot(positions, pinv(slab.cell))

        # Convert to fractional
        unique_positions, unique = [], []
        for i, xyz in enumerate(frac_coords):

            symmetry_match = False
            for j, rotation in enumerate(rotations):
                translation = translations[j]

                affine_matrix = np.eye(4)
                affine_matrix[0:3][:, 0:3] = rotation
                affine_matrix[0:3][:, 3] = translation

                affine_point = np.array([xyz[0], xyz[1], xyz[2], 1])
                operation = np.dot(affine_matrix, affine_point)[0:3]

                if len(matching_sites(operation, unique_positions)) > 0:
                    symmetry_match = True
                    break

            if not symmetry_match:
                unique_positions += [xyz]
                unique += [i]

        unique = np.array(unique)
        sites[k][0] = np.dot(unique_positions, slab.cell)
        sites[k][1] = np.array(points[unique])

        if isinstance(nnneighbors, np.ndarray):
            sites[k][2] = nnneighbors[unique]

    return sites


def matching_sites(position, comparators, tol=1e-8):
    """ Get the indices of all points in a comparator list that are
    equal to a fractional coordinate (with a tolerance), taking into
    account periodic boundary conditions (adaptation from Pymatgen).

    Parameters:
      position: list (3,)
        Fractional coordinate to compare to list.
      comparators: list (3, n)
        Fractional coordinates to compare against.
      tol: float
        Absolute tolerance.

    Returns: list (n,)
        Indices of matches.
    """

    if len(comparators) == 0:
        return []

    fdist = comparators - position
    fdist -= np.round(fdist)
    return np.where((np.abs(fdist) < tol).all(axis=1))[0]


def _get_adsorption_vectors(
        vslab,
        sites,
):
    """ Returns the vectors representing the furthest distance from
    the neighboring atoms.

    (TODO: This input is complex and confusing. Would be nice to simplify.)

    Parameters:
      vslab: ASE atoms-object
        The virtual surface produced from find_adsorption_sites.
      sites: dict of 3 lists
        Dictionary of top, bridge, hollow, and 4-fold sites containing
        positions, points, and neighbor lists.

    Returns: dict of 3 lists
      Dictionary of top, bridge, hollow, and 4-fold sites containing
      positions, points, and adsorption vector lists.
    """

    pos = vslab.positions
    for k, v in sites.items():
        coordinates, points, neighbors = v

        vectors = []
        for i, s in enumerate(coordinates):

            if len(neighbors):
                direct = pos[np.append(points[i], neighbors[i])]
            else:
                direct = pos[points[i]]

            vectors += [utils.plane_normal(direct)]

        sites[k][2] = np.array(vectors)

    return sites
