from scipy.spatial import Voronoi
from ase.data import covalent_radii as radii
import numpy as np
import spglib
from ase import Atoms


def expand_cell(atoms, r=6):
    """ Return Cartesian coordinates atoms within a supercell
    which contains spheres of specified cutoff radius around
    all atom positions.

    Parameters:
      atoms: ASE atoms-object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.

      r: float
        Radius of the spheres to expand around each atom.

    Returns:
      index: ndarray of int
        Indices associated with the original unit cell positions.

      coords: ndarray of (3,) array
        Cartesian coordinates associated with positions in the
        supercell.
    """

    cell = atoms.get_cell()
    recp_len = np.diag(np.linalg.pinv(cell))
    nmax = float(r) * recp_len + 0.01

    pbc = atoms.get_pbc()
    low = np.floor(-nmax * pbc)
    high = np.ceil(nmax * pbc + 1)

    arange = np.arange(low[0], high[0])
    brange = np.arange(low[1], high[1])
    crange = np.arange(low[2], high[2])

    arange = arange[:, None] * np.array([1, 0, 0])[None, :]
    brange = brange[:, None] * np.array([0, 1, 0])[None, :]
    crange = crange[:, None] * np.array([0, 0, 1])[None, :]

    images = arange[:, None, None] + \
             brange[None, :, None] + \
             crange[None, None, :]

    cart_images = np.dot(images, cell)

    coords = atoms.positions[:, None, None, None, :] + \
             cart_images[None, :, :, :, :]

    index = np.where(coords.sum(axis=4))[0]
    index = np.append([0], index)

    coords = coords.flatten()
    n = int(coords.shape[0] / 3)
    coords = coords.reshape((n, 3))

    return index, coords


def get_voronoi_neighbors(atoms, r=10):
    """ Return the nearest-neighbors list from the Voronoi
    method.

    Parameters:
      atoms: ASE atoms-object
        Atoms object with the periodic boundary conditions and
        unit cell information to use.

      r: float
        Radius of the spheres to expand around each atom.

    Returns:
      neighbors: dict
        Edge tuples notation denoting bonds between the atoms of
        the corresponding indices and values of the number of
        those bonds.

        Multi-bonding occurs through periodic boundary conditions.
    """

    index, coords = expand_cell(atoms, r)

    oind = np.empty(len(atoms))
    for i, A in enumerate(atoms.positions):
        for j, B in enumerate(coords):
            if np.allclose(A, B):
                oind[i] = j
                break

    voronoi = Voronoi(coords)
    edges = {}
    indices = np.zeros(len(atoms))
    for nn, vind in voronoi.ridge_dict.items():
        if -1 in vind:
            continue

        for n in nn:
            if n in oind:
                indices[index[n]] += 1
                e = tuple(sorted((index[nn[0]], index[nn[1]])))
                if e not in edges:
                    edges[e] = 0.5
                else:
                    edges[e] += 0.5

    return indices, edges


def get_cutoff_neighbors(atoms, cutoff=4, atol=1e-8):
    """ I can haz documentation?
    """

    cutoff = cutoff + atol

    index, coords = expand_cell(atoms, cutoff * 2)

    indices = np.zeros(len(atoms))
    edges = {}
    for i, center in enumerate(atoms.positions):
        norm = coords - center
        d2 = (norm ** 2).sum(axis=1)
        srt = np.argsort(d2)

        d2 = d2[srt][1:]
        nindex = index[srt][1:]

        inner = np.where(d2 <= cutoff ** 2)[0]

        for j in nindex[inner]:
            indices[i] += 1
            e = tuple(sorted((i, j)))
            if e not in edges:
                edges[e] = 0.5
            else:
                edges[e] += 0.5

    return indices, edges


def get_neighbors(
        atoms,
        points=None,
        cutoff_matrix=None):
    """ Returns the neighboring atoms within a specified cutoff matrix
    criteria for requested points.

    Use of the cutoff matrix provides more fine-tuned control
    over the interaction parameters.

    Parameters:
      atoms: ASE atoms-object
        Atoms object to return.

      points: list (N,) or None
        Points to locate neighboring points of. If not provided, all
        atom points in the atoms-object will be used.

      cutoff_matrix: ndarray (96, 96) or None
        A matrix of interaction distances for the first 96 elements
        of the periodic table. These interactions are separated into
        individual i, j interactions. If None, defaults from ASE
        will be used.

    Returns:
      neighbors: dict
        Keys of each point and an array of each neighboring atoms index.
    """

    if cutoff_matrix is None:

        r = radii.copy()
        # TODO: develop an SVM to parameterize this for me
        # Will need reliable training data or an unsupervised approach
        metals = [
            [47, 1.1],  # Ag
            [79, 1.2],  # Au
            [29, 1.1],  # Cu
            [77, 1.1],  # Ir
            [46, 1.1],  # Pd
            [78, 1.2],  # Pt
            [45, 1.1],  # Rh
        ]

        adsorbates = [
            [1, 1.0],
            [6, 1.0],
            [7, 1.0],
            [8, 1.0],
            [16, 1.0]
        ]

        for i, f in metals:
            r[i] *= f
        for i, f in adsorbates:
            r[i] *= f
        cutoff_matrix = np.zeros((96, 96))
        for i in range(96):
            for j in range(96):
                cutoff_matrix[i][j] = r[i] + r[j]
                cutoff_matrix[j][i] = r[i] + r[j]

    rcmax = cutoff_matrix.max()

    an = atoms.get_atomic_numbers()
    positions = atoms.get_positions()
    pbc = atoms.get_pbc()
    cell = atoms.get_cell()

    icell = np.linalg.pinv(cell)
    scaled = np.dot(positions, icell)
    scaled0 = scaled.copy()

    N = []
    for i in range(3):
        if pbc[i]:
            scaled0[:, i] %= 1.0
            v = icell[:, i]
            h = 1 / np.sqrt(np.dot(v, v))
            n = int(2 * rcmax / h) + 1
        else:
            n = 0
        N.append(n)

    offsets = (scaled0 - scaled).round().astype(int)
    positions0 = atoms.positions + np.dot(offsets, cell)
    natoms = len(atoms)
    indices = np.arange(natoms)

    if points is None:
        points = indices

    cutoffs = np.zeros(natoms)
    neighbors = {a: np.empty(0, int) for a in points}
    for n1 in range(-N[0], N[0] + 1):
        for n2 in range(-N[1], N[1] + 1):
            for n3 in range(-N[2], N[2] + 1):

                displacement = np.dot((n1, n2, n3), cell)

                for a in points:

                    for b in range(natoms):
                        cutoffs[b] = cutoff_matrix[an[a]][an[b]]

                    d = positions0 + displacement - positions0[a]
                    i = indices[(d**2).sum(1) < (cutoffs)**2]

                    if a in i:
                        i = np.delete(i, np.where(i == a)[0])

                    neighbors[a] = np.concatenate(
                        (neighbors[a], i))

    return neighbors


def get_primitive_cell(
        atoms,
        tol=1e-8):
    """ ASE atoms-object interface with spglib primitive cell finder:
    https://atztogo.github.io/spglib/python-spglib.html#python-spglib

    Parameters:
      atoms: ASE atoms-object
        Atoms object to search for a primitive unit cell.

      tol: float
        Tolerance for floating point rounding errors.

    Returns:
      primitive cell: ASE atoms-object
        The primitive unit cell returned by spglib if one is found.
    """

    lattice = atoms.cell
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()

    cell = (lattice, positions, numbers)

    _lattice, _positions, _numbers = spglib.find_primitive(
        cell,
        symprec=tol)

    atoms = Atoms(symbols=_numbers, cell=_lattice)
    atoms.set_scaled_positions(_positions)

    return atoms


def get_symmetry(
        atoms,
        tol=1e-8):
    """ ASE atoms-object interface with spglib symmetry finder:
    https://atztogo.github.io/spglib/python-spglib.html#python-spglib

    Parameters:
      atoms: ASE atoms-object
        Atoms object to search for symmetric structures of.

      tol: float
        Tolerance for floating point rounding errors.

    Returns:
      symmetry operations: ndarray (N, N)
        Symmetry operations from spglib.
    """

    lattice = atoms.cell
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()

    cell = (lattice, positions, numbers)

    return spglib.get_symmetry(cell, symprec=tol)


def get_unique_coordinates(
        atoms,
        axis=2,
        direct=True,
        tag=False,
        tol=1e-5):
    """ Return unique coordinate values of a given atoms object
    for a specified axis.

    Parameters:
      atoms: ASE atoms-object
        Atoms object to search for unique values along.

      axis: int
        Value of 0, 1, or 2 associated with x, y, and z coordinates.

      direct: bool
        Whether to use direct coordinates or Cartesian.

      tag: bool
        Assign ase-like tags to each layer of the slab.

      tol: float
        The tolerance to search for unique values within.

    Returns:
      values: ndarray (N,)
        Array of unique values.
    """

    if direct:
        positions = atoms.get_scaled_positions()
    else:
        positions = atoms.positions

    values = [positions[0][axis]]
    for d in positions[1:, axis]:
        if not np.isclose(d, values, rtol=tol).any():
            values += [d]
    values = np.sort(values)

    if tag:
        tags = []
        for p in positions[:, axis]:
            close = np.isclose(p, values[::-1], rtol=tol)
            tags += [np.where(close)[0][0] + 1]
        atoms.set_tags(tags)

    return values


# tmp for adsorption testing
from ase.io import read
from ase import Atoms, Atom
from ase.constraints import FixAtoms
from ase.data import covalent_radii, atomic_numbers


def read_structure(filename):
    '''
      Description: Read the alat and basis from given path...
      Input:
        filename: traj file to be read
      Output:
        alat: lattice vectors
        basis: basis for structure
    '''

    alat = None
    basis = None

    out = read(filename, index=-1)
    symbols = out.get_chemical_symbols()
    positions = out.get_positions()
    cell = out.get_cell()
    _tags = []

    try:
        _tags = out.constraints[0].index
    except BaseException:
        pass

    tags = [-1 for s in positions]
    for i in _tags:
        tags[i] = 1

    alat = cell
    inv_alat = np.linalg.inv(alat)
    basis = []
    for i in range(len(positions)):
        symbol = symbols[i]
        pos_cart = positions[i]
        pos_direct = np.matmul(pos_cart, inv_alat)
        if_fix = tags[i]
        basis.append({'symbol': symbol, 'pos_cart': np.array(
            pos_cart), 'pos_direct': np.array(pos_direct), 'if_fix': if_fix})

    return alat, basis


def get_closest_image_of_atom(atom_i, atom_j, alat):
    ''' Return the minimum distance image (by considering PBC)
    of atom_j by considering atom_i as center.
    '''

    dist = 1e+6
    image = [0, 0, 0]
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                pos_direct = np.add(atom_j['pos_direct'], [i, j, k])
                pos_cart = np.matmul(pos_direct, alat)
                dist_ = np.sum(
                    (np.add(atom_i['pos_cart'], -pos_cart))**2)
                if(dist_ < dist):
                    dist = dist_
                    image = [i, j, k]

    pos_direct = np.add(atom_j['pos_direct'], image)
    pos_cart = np.matmul(pos_direct, alat)
    atom = {
        'symbol': atom_j['symbol'],
        'if_fix': atom_j['if_fix'],
        'pos_direct': pos_direct,
        'pos_cart': pos_cart}

    return atom


def get_distance_between_atoms(atom_i, atom_j):
    "return the distance between two atoms"

    dist = np.sqrt(np.sum(
        np.add(atom_i['pos_cart'], -np.array(atom_j['pos_cart']))**2))

    return dist


def get_cartesian_coordinate_from_direct_coordinate(alat, _basis):
    '''
      Description: given basis, populate cartesian positions from direct
      Input:
        alat: lattice_vectors
        _basis: basis to convert
      output:
        basis: converted basis
    '''

    basis = []
    for b in _basis:
        pos_direct = b['pos_direct']
        pos_cart = np.matmul(pos_direct, alat)
        basis.append({'symbol': b['symbol'],
                      'pos_cart': np.array(pos_cart),
                      'pos_direct': np.array(pos_direct),
                      'if_fix': b['if_fix']})

    return basis


def get_center_of_structure(alat, basis):
    '''
      Description: find the cneter of given structure
      Input:
        alat: lattice vectors of structure
        basis: baiss of structure
      Output:
        [x,y,z] cartesian coordinate of cenetr given structure
    '''

    center = [0, 0, 0]
    for b in basis:
        center = np.add(center, b['pos_cart'])
    center = center * (1.0 / len(basis))

    return center


def get_bottom_of_structure(alat, basis):
    '''
      Description: find the bottom of given structure
      Input:
        alat: lattice vectors of structure
        basis: baiss of structure
      Output:
        z-coordinate (in cartesian) of bottom of given structure
    '''

    bottom = 1e+6
    for b in basis:
        bottom = min(bottom, b['pos_cart'][2])

    return bottom


def get_direct_coordinate_from_cartesian_coordinate(alat, _basis):
    '''
      Description: given basis, populate direct positions from cart positions
      Input:
        alat: lattice_vectors
        _basis: basis to convert
      output:
        basis: converted basis
    '''

    inv_alat = np.linalg.inv(alat)
    basis = []
    for b in _basis:
        pos_cart = b['pos_cart']
        pos_direct = np.matmul(pos_cart, inv_alat)
        basis.append({'symbol': b['symbol'],
                      'pos_cart': np.array(pos_cart),
                      'pos_direct': np.array(pos_direct),
                      'if_fix': b['if_fix']})

    return basis


def get_ase_atoms_object(alat, basis):
    # generae ase atoms object

    atoms = Atoms([Atom(atom['symbol'],
                        (atom['pos_cart']),
                        tag=atom['if_fix']) for atom in basis])
    atoms.set_cell(np.array(alat))
    atoms.set_pbc(['True', 'True', 'True'])
    constraint = FixAtoms(
        indices=[atom.index for atom in atoms if atom.tag == 1])
    atoms.set_constraint(constraint)

    return atoms


def get_atomic_radius(atom):
    "return the radius of atom"

    return covalent_radii[atomic_numbers[atom]]


def get_surface_supercell(alat, basis):
    "return the surface supercell with size [3*3*1]"

    supercell = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            for b in basis:
                supercell.append({
                    'symbol': b['symbol'],
                    'pos_direct': np.add(b['pos_direct'], [i, j, 0]),
                    'size': get_atomic_radius(b['symbol']),
                    'pos_cart': [],
                    'if_fix': b['if_fix']})

    for b in supercell:
        b['pos_cart'] = np.matmul(b['pos_direct'], alat)

    return supercell


def get_rotation_matrix(u, theta):
    "return rotation matrix around u with angle theta. theta is in radians"

    ct = float(np.cos(theta))
    st = float(np.sin(theta))
    ux = u[0]
    uy = u[1]
    uz = u[2]

    R = [[ct + ux * ux * (1 - ct),
          ux * uy * (1 - ct) - uz * st,
          ux * uz * (1 - ct) + uy * st],
         [uy * ux * (1 - ct) + uz * st,
          ct + uy * uy * (1 - ct),
          uy * uz * (1 - ct) - ux * st],
         [uz * ux * (1 - ct) - uy * st,
          uz * uy * (1 - ct) + ux * st,
          ct + uz * uz * (1 - ct)]]

    return R


def sortkey_distance(item):
    "sort the list by dist"

    return item['dist']


def if_same_atom_based_on_neighbors(neighbor_i, neighbor_j, symprec=1e-8):
    '''
      Description: check if given neighborlists corresponds to same atom
      Input:
        neighbors_i, neighbors_j: neighbor_list of atom_i, atom_j
        symprec: precision for floating point math
      Output:
        True if same neighbor_list
    '''

    # print(len(neighbor_i))

    if(len(neighbor_i) != len(neighbor_j)):
        return False

    # check if all neighbor of i are in j
    for ii in range(len(neighbor_i)):
        if_neighbor_found = 0

        for jj in range(len(neighbor_j)):

            if(neighbor_i[ii]['symbol'] != neighbor_j[jj]['symbol']):
                continue
            if(abs(neighbor_i[ii]['dist'] - neighbor_j[jj]['dist']) > symprec):
                continue
            if_neighbor_found = 1
            break

        if(if_neighbor_found == 0):
            return False

    # check if all neighbors of j are in i
    for jj in range(len(neighbor_j)):
        if_neighbor_found = 0

        for ii in range(len(neighbor_i)):

            if(neighbor_i[ii]['symbol'] != neighbor_j[jj]['symbol']):
                continue
            if(abs(neighbor_i[ii]['dist'] - neighbor_j[jj]['dist']) > symprec):
                continue
            if_neighbor_found = 1
            break

        if(if_neighbor_found == 0):
            return False

    return True


def get_neighbors_of_given_atom(
        atom,
        basis,
        neighbor_cutoff=4.0,
        symprec=1e-8):
    ''' find all neighbors of given atom in the given basis within
    neighbor_cutoff distance'''

    neighbor = []
    for _atom in basis:
        dist = np.sum(np.add(atom['pos_cart'], -_atom['pos_cart'])**2)

        if(dist < (neighbor_cutoff + symprec)**2):
            neighbor.append({'symbol': _atom['symbol'], 'dist': dist, 'pos_direct': np.array(
                _atom['pos_direct']), 'pos_cart': np.array(_atom['pos_cart'])})

    neighbor.sort(key=sortkey_distance)
    return neighbor
