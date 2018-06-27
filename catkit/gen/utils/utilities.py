from catkit import Gratoms
from ase.data import chemical_symbols as sym
import numpy as np
import re
try:
    from math import gcd
except ImportError:
    from fractions import gcd


def running_mean(array, N=5):
    """Calculate the running mean of array for N instances.

    Parameters
    ----------
    array : array_like | ndarray (N,)
        Array of values to have a average taken from.
    N : int
        Number of values to take an average with.

    Returns
    -------
    running_mean : ndarray (N + 1,)
        Mean value of the running average.
    """
    length = len(array)
    if length < N:
        N = length

    cumsum = np.cumsum(np.insert(array, 0, 0))
    running_mean = (cumsum[N:] - cumsum[:-N]) / float(N)

    return running_mean


def to_gratoms(atoms):
    """Convert and atom object to a gratoms object."""
    gratoms = Gratoms(
        numbers=atoms.numbers,
        positions=atoms.positions,
        pbc=atoms.pbc,
        cell=atoms.cell)

    if atoms.constraints:
        gratoms.set_constraint(atoms.constraints)

    return gratoms


def get_atomic_numbers(formula, return_count=False):
    """Return the atomic numbers associated with a chemical formula.

    Parameters
    ----------
    formula : string
        A chemical formula to parse into atomic numbers.
    return_count : bool
        Return the count of each element in the formula.

    Returns
    -------
    numbers : ndarray (n,)
        Element numbers in associated species.
    counts : ndarray (n,)
        Count of each element in a species.
    """
    parse = re.findall('[A-Z][a-z]?|[0-9]+', formula)

    values = {}
    for i, e in enumerate(parse):
        if e.isdigit():
            values[parse[i - 1]] += int(e) - 1
        else:
            if e not in values:
                values[e] = 1
            else:
                values[e] += 1

    numbers = np.array([sym.index(k) for k in values.keys()])
    srt = np.argsort(numbers)
    numbers = numbers[srt]

    if return_count:
        counts = np.array([v for v in values.values()])[srt]

        return numbers, counts

    return numbers


def get_reference_energies(species, energies):
    """Get reference energies for the elements in a set of molecules.

    Parameters
    ----------
    species : list (n,)
        Chemical formulas for each molecular species.
    energies : list (n,)
        Total energies associated with each species.

    Returns
    -------
    elements : ndarray (n,)
        Atomic elements associated with all species.
    references : ndarray (n,)
        Reference energies associated with each element.
    """
    if not isinstance(energies, np.ndarray):
        energies = np.array(energies)

    A = np.zeros((len(species), len(species)))
    elements = np.zeros(len(species), dtype=int)
    n = 0

    # Construct the elements array as they appear
    for i, s in enumerate(species):
        num, cnt = get_atomic_numbers(s, True)

        for j in num[~np.in1d(num, elements)]:
            elements[n] = j
            n += 1

        A[i][np.in1d(elements, num)] = cnt

    references = np.linalg.solve(A, energies)
    srt = np.argsort(elements)
    references = references[srt]
    elements = elements[srt]

    return elements, references


def parse_slice(slice_name):
    """Return a correctly parsed slice from input of varying types."""
    if isinstance(slice_name, (slice)):
        _slice = slice_name

    elif isinstance(slice_name, type(None)):
        _slice = slice(None)

    elif isinstance(slice_name, int):
        i = int(slice_name)
        _slice = slice(i, i + 1)

    elif isinstance(slice_name, str):
        if slice_name.isdigit():
            i = int(slice_name)
            _slice = slice(i, i + 1)

        else:
            split = slice_name.split(':')
            split = [int(_) if _.lstrip('-').isdigit()
                     else None for _ in split]
            _slice = slice(*split)

    return _slice


def ext_gcd(a, b):
    """Extension of greatest common divisor."""
    if b == 0:
        return 1, 0
    elif a % b == 0:
        return 0, 1
    else:
        x, y = ext_gcd(b, a % b)
        return y, x - y * (a // b)


def list_gcd(values):
    """Return the greatest common divisor of a list of values."""
    if isinstance(values[0], float):
        values = np.array(values, dtype=int)

    gcd_func = np.frompyfunc(gcd, 2, 1)
    _gcd = np.ufunc.reduce(gcd_func, values)

    return _gcd
