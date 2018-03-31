import numpy as np
from numpy.linalg import matrix_rank, det
from itertools import combinations
import math


def get_independent_set(nu):
    """Returns an array of linearly independent reaction routes.

    Parameters:
    -----------
    nu : ndarray (n, m)
        The stoichiometric matrix of n species and m mechanisms.

    Returns:
    --------
    sigma : ndarray (l, n)
        A linearly independent set of possible reaction routes.
    """
    p = nu.shape[0]
    m = matrix_rank(nu)
    reactions = np.arange(p)
    selection = reactions == 0

    for n in reactions[~selection]:
        selection[n] = True
        rank = matrix_rank(nu[selection, :m])
        if rank == len(reactions[selection]):
            continue

        selection[n] = False

    # Find a linearly independent set of RR
    sigma = np.zeros((p - m, p), dtype=int)
    for i, n in enumerate(reactions[~selection]):
        selection[n] = True

        _nu = nu[selection, :m][None]
        S = np.repeat(_nu, len(reactions[selection]), axis=0)
        eye = np.eye(len(reactions[selection]))[:, :, None]
        R = np.concatenate([S, eye], axis=2)

        # Does not convert to correct integers without round
        values = np.round(det(R)).astype(int)

        # Screen trivial solutions
        if np.all(values == 0):
            selection[n] = False
            continue

        # Normalize first index as positive and reduce by gcd
        _gcd = list_gcd(values)
        values *= np.sign(values[np.nonzero(values)[0][0]])
        values = (values / _gcd).astype(int)

        sigma[i, selection] = values
        selection[n] = False

    return sigma


def get_response_reactions(epsilon, selection=None, species=False):
    """Returns an array of possible response reaction for a given
    chemical formula array.

    Parameters:
    -----------
    epsilon : ndarray (n, m)
        The chemical formula array of n elements by m molecular species.
    species : bool
        Return the indices of the chemical species used

    Returns:
    --------
    RER : ndarray (l, m)
        Possible response reactions.
    index : ndarray (l, k)
        Indices of the k chemical species use to produce the
        l response reactions.
    """
    s = matrix_rank(epsilon)
    RER, index = [], []
    if not selection:
        selection = np.arange(epsilon.shape[0])

    for sel in combinations(selection, r=s + 1):
        values = np.zeros(epsilon.shape[0], dtype=int)

        sigma = np.repeat(epsilon[[sel]][None], s + 1, axis=0)
        eye = np.eye(s + 1)[:, :, None]
        R = np.concatenate([sigma, eye], axis=2)

        # Does not convert to correct integers without round
        values[[sel]] = np.round(det(R))

        # Screen trivial solutions
        if np.all(values == 0):
            continue

        # Normalize first index as positive and reduce by gcd
        _gcd = list_gcd(values)
        values *= np.sign(values[np.nonzero(values)[0][0]])
        values = (values / _gcd).astype(int)

        # Screen the stoichiometric matches
        match = False
        for v in RER:
            match = np.all(v == values)
            if match:
                break

        if not match:
            RER += [values]
            index += [list(sel)]

    RER = np.array(RER)
    if species:
        index = np.array(index)
        return RER, index

    return RER


def list_gcd(values):
    """Return the greatest common denominator of values in a list."""
    gcd = np.frompyfunc(math.gcd, 2, 1)
    list_gcd = np.ufunc.reduce(gcd, values)

    return list_gcd
