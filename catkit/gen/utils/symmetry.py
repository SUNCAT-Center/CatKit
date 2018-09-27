from ..classify import get_modified_spin_symbols
from ..classify import get_unmodified_spin_symbols
from catkit import Gratoms
import numpy as np
import spglib


def get_spglib_cell(atoms, primitive=False, tol=1e-5):
    """Atoms object interface with spglib primitive cell finder:
    https://atztogo.github.io/spglib/python-spglib.html#python-spglib

    The function also builds in limited functionality for initial
    magnetic moments. Only integer values are supported.

    Parameters
    ----------
    atoms : object
        Atoms object to search for a primitive unit cell.
    primitive : bool
        Reduce the atoms object into a primitive form.
    tol : float
        Tolerance for floating point rounding errors.

    Returns
    -------
    primitive cell : object
        The primitive unit cell returned by spglib if one is found.
    """
    lattice = atoms.cell
    positions = atoms.get_scaled_positions()
    numbers = atoms.numbers
    magmoms = atoms.get_initial_magnetic_moments()
    modified_numbers = get_modified_spin_symbols(numbers, magmoms)

    cell = (lattice, positions, modified_numbers)
    cell = spglib.standardize_cell(cell, to_primitive=primitive, symprec=tol)

    if cell is None:
        return atoms

    _lattice, _positions, _modified_numbers = cell
    _numbers, _magmoms = get_unmodified_spin_symbols(_modified_numbers)

    atoms = Gratoms(symbols=_numbers, cell=_lattice, pbc=atoms.pbc)
    atoms.set_scaled_positions(_positions)
    atoms.set_initial_magnetic_moments(_magmoms)

    return atoms


class Symmetry():
    """Wrapper for the spglib package."""

    def __init__(self, atoms, tol=1e-5):
        """Atoms object interface with spglib symmetry finder:
        https://atztogo.github.io/spglib/python-spglib.html#python-spglib

        Parameters
        ----------
        atoms : Atoms object
            Atomic structure to return the symmetry operations for.
        tol : float
            Tolerance for floating point rounding errors.
        """
        self.lattice = atoms.cell
        self.positions = atoms.get_scaled_positions()
        self.numbers = atoms.get_atomic_numbers()
        self.magmoms = atoms.get_initial_magnetic_moments()
        self.modified_numbers = get_modified_spin_symbols(
            self.numbers, self.magmoms)
        self.tol = tol

        cell = (self.lattice, self.positions, self.modified_numbers)
        self.data = spglib.get_symmetry_dataset(cell, symprec=tol)

    def get_symmetry(self, operations=True):
        """Return the symmetry operations for a given atomic structure.

        Parameters
        ----------
        operations : bool
            Whether to return the affine matrix operations.

        Returns
        -------
        rotations : ndarray (N, 3, 3)
            Rotation matices of the symmetry operations.
        translations ndarray (N, 3)
            Translation vector components of the symmetry operations.
        affine_matrices ndarray (N, 4, 4)
            Affine matrix operations, combinations of the rotation and
            translation with ones along the diagonal.
        """
        rotations = self.data['rotations']
        translations = self.data['translations']

        if operations:
            affine_matrices = np.zeros((rotations.shape[0], 4, 4))
            affine_matrices[:, :3, :3] = rotations
            affine_matrices[:, -1, :3] = translations
            affine_matrices[:, -1, -1] = 1
            return affine_matrices

        return rotations, translations

    def get_pointgroup(self, check_laue=False):
        """Return the point group operations of a systems.

        Parameters
        ----------
        check_laue : bool
            Return if the pointgroup is a laue symmetry.

        Returns
        -------
        pointgroup : str
            The pointgroup symmetry of the atomic structure.
        is_laue : bool
            Whether the pointgroup is a laue symmetry.
        """
        pointgroup = self.data['pointgroup']

        if check_laue:
            laue = ['-1', '2/m', 'mmm', '4/m', '4/mmm',
                    '-3', '-3m', '6/m', '6/mmm', 'm-3', 'm-3m']
            is_laue = pointgroup in laue

            return pointgroup, is_laue

        return pointgroup

    def get_lattice_name(self):
        """Return the lattice name of an atoms object based
        on its spacegroup number:
        https://en.wikipedia.org/wiki/List_of_space_groups

        Returns
        -------
        lattice : str
            The name of the structures lattice.
        """
        space_group_number = self.data['number']

        if space_group_number in [146, 148, 155, 160, 161, 166, 167]:
            return 'rhombohedral'

        lattices = {
            'triclinic': 2,
            'monoclinic': 15,
            'orthorhombic': 74,
            'tetragonal': 142,
            'hexagonal': 194,
            'cubic': 230}

        for lattice, max_number in lattices.items():
            if space_group_number <= max_number:
                return lattice

    def get_conventional_standard_cell(self):
        """This function is copied from Pymatgen.

        Gives a structure with a conventional cell according to certain
        standards. The standards are defined in Setyawan, W., & Curtarolo,
        S. (2010). High-throughput electronic band structure calculations:
        Challenges and tools. Computational Materials Science,
        49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010
        They basically enforce as much as possible
        norm(a1)<norm(a2)<norm(a3)

        Parameters
        ----------
        atoms : Atoms object
            The atoms object to return the conventional standard form of.
        tol : float
            Float point precision tolerance.

        Returns
        -------
            The structure in a conventional standardized cell.
        """
        abc = np.linalg.norm(self.lattice, axis=1)
        lattice_name = self.get_lattice_name()
        spacegroup = self.data['international']
        transform = np.zeros((3, 3))

        if lattice_name in ['cubic', 'orthorhombic']:
            if spacegroup.startswith('C'):
                sort = np.argsort(abc[:2])
                abc[[0, 1], :] = abc[sort, :]
                transform[[0, 1], sort] = 1
                transform[2][2] = 1

            elif spacegroup.startswith('A'):
                sort = np.argsort(abc[1:])
                c = abc[0].copy()
                abc[[1, 2], :] = abc[sort, :]
                abc[0] = c
                transform[[0, 1], sort] = 1
                transform[2][0] = 1
            else:
                sort = np.argsort(abc)
                transform[[0, 1, 2], sort] = 1
                abc = abc[sort]
            # latt = Lattice.orthorhombic(a, b, c)

        elif lattice_name == 'tetragonal':
            sort = np.argsort(abc)
            transform[[0, 1, 2], sort] = 1

            if abs(abc[1] - abc[2]) < self.tol and \
               abs(abc[0] - abc[2]) > self.tol:
                abc[[0, 2], :] = abc[[2, 0], :]
                transform = np.dot([[0, 0, 1], [0, 1, 0], [1, 0, 0]],
                                   transform)
            latt = Lattice.tetragonal(a, c)

        elif lattice_name in ['hexagonal', 'rhombohedral']:
            # check first if we have the refined structure shows a rhombohedral
            # cell
            # if so, make a supercell
            a, b, c = latt.abc
            if np.all(np.abs([a - b, c - b, a - c]) < 0.001):
                struct.make_supercell(((1, -1, 0), (0, 1, -1), (1, 1, 1)))
                a, b, c = sorted(struct.lattice.abc)

            if abs(b - c) < 0.001:
                a, c = c, a
            new_matrix = [[a / 2, -a * math.sqrt(3) / 2, 0],
                          [a / 2, a * math.sqrt(3) / 2, 0],
                          [0, 0, c]]
            latt = Lattice(new_matrix)
            transf = np.eye(3, 3)

        elif latt_type == 'monoclinic':
            # You want to keep the c axis where it is to keep the C- settings

            if self.get_space_group_operations().int_symbol.startswith('C'):
                transf = np.zeros(shape=(3, 3))
                transf[2] = [0, 0, 1]
                sorted_dic = sorted([{'vec': latt.matrix[i],
                                      'length': latt.abc[i],
                                      'orig_index': i} for i in [0, 1]],
                                    key=lambda k: k['length'])
                a = sorted_dic[0]['length']
                b = sorted_dic[1]['length']
                c = latt.abc[2]
                new_matrix = None
                for t in itertools.permutations(list(range(2)), 2):
                    m = latt.matrix
                    landang = Lattice(
                        [m[t[0]], m[t[1]], m[2]]).lengths_and_angles
                    if landang[1][0] > 90:
                        # if the angle is > 90 we invert a and b to get
                        # an angle < 90
                        landang = Lattice(
                            [-m[t[0]], -m[t[1]], m[2]]).lengths_and_angles
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = -1
                        transf[1][t[1]] = -1
                        transf[2][2] = 1
                        a, b, c = landang[0]
                        alpha = math.pi * landang[1][0] / 180
                        new_matrix = [[a, 0, 0],
                                      [0, b, 0],
                                      [0, c * cos(alpha), c * sin(alpha)]]
                        continue

                    elif landang[1][0] < 90:
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = 1
                        transf[1][t[1]] = 1
                        transf[2][2] = 1
                        a, b, c = landang[0]
                        alpha = math.pi * landang[1][0] / 180
                        new_matrix = [[a, 0, 0],
                                      [0, b, 0],
                                      [0, c * cos(alpha), c * sin(alpha)]]

                if new_matrix is None:
                    # this if is to treat the case
                    # where alpha==90 (but we still have a monoclinic sg
                    new_matrix = [[a, 0, 0],
                                  [0, b, 0],
                                  [0, 0, c]]
                    transf = np.zeros(shape=(3, 3))
                    for c in range(len(sorted_dic)):
                        transf[c][sorted_dic[c]['orig_index']] = 1
            #if not C-setting
            else:
                # try all permutations of the axis
                # keep the ones with the non-90 angle=alpha
                # and b<c
                new_matrix = None
                for t in itertools.permutations(list(range(3)), 3):
                    m = latt.matrix
                    landang = Lattice(
                        [m[t[0]], m[t[1]], m[t[2]]]).lengths_and_angles
                    if landang[1][0] > 90 and landang[0][1] < landang[0][2]:
                        landang = Lattice(
                            [-m[t[0]], -m[t[1]], m[t[2]]]).lengths_and_angles
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = -1
                        transf[1][t[1]] = -1
                        transf[2][t[2]] = 1
                        a, b, c = landang[0]
                        alpha = math.pi * landang[1][0] / 180
                        new_matrix = [[a, 0, 0],
                                      [0, b, 0],
                                      [0, c * cos(alpha), c * sin(alpha)]]
                        continue
                    elif landang[1][0] < 90 and landang[0][1] < landang[0][2]:
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = 1
                        transf[1][t[1]] = 1
                        transf[2][t[2]] = 1
                        a, b, c = landang[0]
                        alpha = math.pi * landang[1][0] / 180
                        new_matrix = [[a, 0, 0],
                                      [0, b, 0],
                                      [0, c * cos(alpha), c * sin(alpha)]]
                if new_matrix is None:
                    # this if is to treat the case
                    # where alpha==90 (but we still have a monoclinic sg
                    new_matrix = [[sorted_lengths[0], 0, 0],
                                  [0, sorted_lengths[1], 0],
                                  [0, 0, sorted_lengths[2]]]
                    transf = np.zeros(shape=(3, 3))
                    for c in range(len(sorted_dic)):
                        transf[c][sorted_dic[c]['orig_index']] = 1

            op = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
            transf = np.dot(op, transf)
            new_matrix = np.dot(op, new_matrix)
            beta = Lattice(new_matrix).beta
            if beta < 90:
                op = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
                transf = np.dot(op, transf)
                new_matrix = np.dot(op, new_matrix)

            latt = Lattice(new_matrix)

        elif latt_type == 'triclinic':
            #we use a LLL Minkowski-like reduction for the triclinic cells
            struct = struct.get_reduced_structure('LLL')

            a, b, c = latt.lengths_and_angles[0]
            alpha, beta, gamma = [math.pi * i / 180
                                  for i in latt.lengths_and_angles[1]]
            new_matrix = None
            test_matrix = [[a, 0, 0],
                          [b * cos(gamma), b * sin(gamma), 0.0],
                          [c * cos(beta),
                           c * (cos(alpha) - cos(beta) * cos(gamma)) /
                           sin(gamma),
                           c * math.sqrt(sin(gamma) ** 2 - cos(alpha) ** 2
                                         - cos(beta) ** 2
                                         + 2 * cos(alpha) * cos(beta)
                                         * cos(gamma)) / sin(gamma)]]

            def is_all_acute_or_obtuse(m):
                recp_angles = np.array(Lattice(m).reciprocal_lattice.angles)
                return np.all(recp_angles <= 90) or np.all(recp_angles > 90)

            if is_all_acute_or_obtuse(test_matrix):
                transf = np.eye(3)
                new_matrix = test_matrix

            test_matrix = [[-a, 0, 0],
                           [b * cos(gamma), b * sin(gamma), 0.0],
                           [-c * cos(beta),
                            -c * (cos(alpha) - cos(beta) * cos(gamma)) /
                            sin(gamma),
                            -c * math.sqrt(sin(gamma) ** 2 - cos(alpha) ** 2
                                           - cos(beta) ** 2
                                           + 2 * cos(alpha) * cos(beta)
                                           * cos(gamma)) / sin(gamma)]]

            if is_all_acute_or_obtuse(test_matrix):
                transf = [[-1, 0, 0],
                          [0, 1, 0],
                          [0, 0, -1]]
                new_matrix = test_matrix

            test_matrix = [[-a, 0, 0],
                           [-b * cos(gamma), -b * sin(gamma), 0.0],
                           [c * cos(beta),
                            c * (cos(alpha) - cos(beta) * cos(gamma)) /
                            sin(gamma),
                            c * math.sqrt(sin(gamma) ** 2 - cos(alpha) ** 2
                                          - cos(beta) ** 2
                                          + 2 * cos(alpha) * cos(beta)
                                          * cos(gamma)) / sin(gamma)]]

            if is_all_acute_or_obtuse(test_matrix):
                transf = [[-1, 0, 0],
                          [0, -1, 0],
                          [0, 0, 1]]
                new_matrix = test_matrix

            test_matrix = [
                [a, 0, 0],
                [-b * cos(gamma), -b * sin(gamma), 0.0],
                [-c * cos(beta),
                 -c * (cos(alpha) - cos(beta) * cos(gamma)) /
                 sin(gamma),
                 -c * math.sqrt(sin(gamma) ** 2 - cos(alpha) ** 2
                                           - cos(beta) ** 2
                                           + 2 * cos(alpha) * cos(beta)
                                           * cos(gamma)) / sin(gamma)]]
            if is_all_acute_or_obtuse(test_matrix):
                transf = [[1, 0, 0],
                          [0, -1, 0],
                          [0, 0, -1]]
                new_matrix = test_matrix

            latt = Lattice(new_matrix)

        new_coords = np.dot(transf, np.transpose(struct.frac_coords)).T
        new_struct = Structure(latt, struct.species_and_occu, new_coords,
                               site_properties=struct.site_properties,
                               to_unit_cell=True)
        return new_struct.get_sorted_structure()


def get_primitive_standard_structure():
    """This function is copied from Pymatgen.

    Gives a structure with a primitive cell according to certain standards
    the standards are defined in Setyawan, W., & Curtarolo, S. (2010).
    High-throughput electronic band structure calculations:
    Challenges and tools. Computational Materials Science,
    49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010

    Returns
    -------
        The structure in a primitive standardized cell
    """
    conv = get_conventional_standard_structure()
    lattice = self.get_lattice_type()

    if 'P' in self.get_space_group_symbol() or lattice == 'hexagonal':
        return conv

    if lattice == 'rhombohedral':
        lengths, angles = conv.lattice.lengths_and_angles
        if abs(lengths[0]-lengths[2]) < 0.0001:
            transf = np.eye
        else:
            transf = np.array([[-1, 1, 1], [2, 1, 1], [-1, -2, 1]],
                              dtype=np.float) / 3

    elif 'I' in self.get_space_group_symbol():
        transf = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]],
                          dtype=np.float) / 2
    elif 'F' in self.get_space_group_symbol():
        transf = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]],
                          dtype=np.float) / 2
    elif 'C' in self.get_space_group_symbol() \
         or 'A' in self.get_space_group_symbol():
        if self.get_crystal_system() == "monoclinic":
            transf = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 2]],
                              dtype=np.float) / 2
        else:
            transf = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 2]],
                              dtype=np.float) / 2
    else:
        transf = np.eye(3)

    new_sites = []
    latt = Lattice(np.dot(transf, conv.lattice.matrix))
    for s in conv:
        new_s = PeriodicSite(
            s.specie, s.coords, latt,
            to_unit_cell=True, coords_are_cartesian=True,
            properties=s.properties)
        if not any(map(new_s.is_periodic_image, new_sites)):
            new_sites.append(new_s)

    if lattice == 'rhombohedral':
        prim = Structure.from_sites(new_sites)
        lengths, angles = prim.lattice.lengths_and_angles
        a = lengths[0]
        alpha = math.pi * angles[0] / 180
        new_matrix = [
            [a * cos(alpha / 2), -a * sin(alpha / 2), 0],
            [a * cos(alpha / 2), a * sin(alpha / 2), 0],
            [a * cos(alpha) / cos(alpha / 2), 0,
             a * math.sqrt(1 - (cos(alpha) ** 2 / (cos(alpha / 2) ** 2)))]]
        new_sites = []
        latt = Lattice(new_matrix)
        for s in prim:
            new_s = PeriodicSite(
                s.specie, s.frac_coords, latt,
                to_unit_cell=True, properties=s.properties)
            if not any(map(new_s.is_periodic_image, new_sites)):
                new_sites.append(new_s)
        return Structure.from_sites(new_sites)

    return Structure.from_sites(new_sites)


def matching_coordinates(position, comparators, tol=1e-8):
    """Get the indices of all points in a comparator list that are
    equal to a given position (with a tolerance), taking into
    account periodic boundary conditions (adaptation from Pymatgen).

    This will only accept a Cartesian coordinate scheme.
    TODO: merge this with matching_sites.

    Parameters
    ----------
    position : list (3,)
        Fractional coordinate to compare to list.
    comparators : list (3, N)
        Fractional coordinates to compare against.
    tol : float
        Absolute tolerance.

    Returns
    -------
    match : list (N,)
        Indices of matches.
    """
    if len(comparators) == 0:
        return []

    fdist = comparators - position[None, :]
    match = np.where((np.abs(fdist) < tol).all(axis=1))[0]

    return match


def get_unique_coordinates(atoms, axis=2, tag=False, tol=1e-3):
    """Return unique coordinate values of a given atoms object
    for a specified axis.

    Parameters
    ----------
    atoms : object
        Atoms object to search for unique values along.
    axis : int (0, 1, or 2)
        Look for unique values along the x, y, or z axis.
    tag : bool
        Assign ASE-like tags to each layer of the slab.
    tol : float
        The tolerance to search for unique values within.

    Returns
    -------
    values : ndarray (n,)
        Array of unique positions in fractional coordinates.
    """
    positions = (atoms.get_scaled_positions()[:, axis] + tol) % 1
    positions -= tol

    values = [positions[0]]
    for d in positions[1:]:
        if not np.isclose(d, values, atol=tol, rtol=tol).any():
            values += [d]
    values = np.sort(values)

    if tag:
        tags = []
        for p in positions:
            close = np.isclose(p, values[::-1], atol=tol, rtol=tol)
            tags += [np.where(close)[0][0] + 1]
        atoms.set_tags(tags)

    return values
