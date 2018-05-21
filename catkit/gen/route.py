import numpy as np
from numpy.linalg import matrix_rank, det
from itertools import combinations
import math


def get_reaction_routes(nu, sigma, empty_routes=True, independent_only=False):
    """Returns an array of reaction routes. Returns all full
    reaction routes by default.

    Parameters
    ----------
    nu : ndarray (n, m)
        The stoichiometric matrix of n species and m mechanisms.
    sigma : ndarray (m, j)
        A linearly independent set of reaction routes.
    empty_routes : bool
        Return the empty routes along with the full routes.
    independent_only : bool
        Return only a linearly independent set of full reaction routes.
        Can take less time.

    Returns
    -------
    FR : ndarray (m, k)
        Enumerated full reaction routes.
    ER : ndarray (m, l)
        Enumerated empty reaction routes.
    """
    p = sigma.shape[1]
    m = sigma.shape[0]

    isFR = ~np.all(np.dot(sigma, nu) == 0, axis=1)
    rsample = np.arange(1, p)

    n = len(np.where(isFR)[0])
    liRR = np.zeros(sigma.shape, dtype=int)
    liRR[:n] = sigma[isFR]

    ER = [rr for rr in sigma[~isFR]]
    FR = [rr for rr in sigma[isFR]]
    for c in combinations(rsample, r=m - 1):
        S = np.repeat(sigma[:, c][None], sigma.shape[0], axis=0)
        eye = np.eye(sigma.shape[0])[:, :, None]
        R = np.concatenate([S, eye], axis=2)

        # Does not convert to correct integers without round
        values = np.round(det(R)).astype(int)

        # Screen trivial solutions
        if np.all(values == 0):
            continue

        # Normalize first index as positive and reduce by gcd
        _gcd = list_gcd(values)
        values = (values / _gcd).astype(int)
        route = (sigma * values[:, None]).sum(axis=0)
        route *= np.sign(route[np.nonzero(route)[0][0]])

        OR = np.dot(route, nu)
        if not np.all(OR == 0):
            match = False
            for rr in FR:
                match = np.all(rr == route)
                if match:
                    break

            if not match:
                FR += [route]

            if independent_only:
                liRR[n] = route
                if matrix_rank(liRR) == n + 1:
                    n += 1

                    if n == m:
                        return liRR
        elif empty_routes:
            match = False
            for rr in ER:
                match = np.all(rr == route)
                if match:
                    break

            if not match:
                ER += [route]

    if empty_routes:
        return np.array(FR), np.array(ER)

    return np.array(FR)


def get_heppel_sellers(nu, terminal):
    """Returns an array of linearly independent reaction routes
    as described by Heppel-Sellers reaction route enumeration.

    Parameters
    ----------
    nu : ndarray (n, m)
        The stoichiometric matrix of n species and m mechanisms.
    terminal : ndarray (j,)
        Indices of the m species to be considered as terminal

    Returns
    -------
    sigma : ndarray (m, k)
        Linearly independent set of Heppel-Sellers reaction routes.
    """
    inter = ~np.in1d(np.arange(nu.shape[1]), terminal)

    # Setting up some commonly used parameters
    n = nu.shape[1]
    m = matrix_rank(nu)
    inter = np.where(inter)[0]
    selection = inter == inter[0]
    nuT = nu.T[inter]

    # Reduce the intermetiates to those which are linearly independent
    for i, s in enumerate(inter[~selection]):
        selection[i + 1] = True

        A = nuT[selection]
        rank = matrix_rank(A)
        if rank == A.shape[0]:
            if rank == m - 1:
                break
            continue

        selection[i + 1] = False

    inter = inter[selection]

    # Now choose m - 1 linearly independent elementary steps
    p = nu.shape[0]
    reactions = np.arange(p)
    selection = reactions == 0

    for n in reactions[~selection]:
        selection[n] = True
        R = nu[selection][:, inter]

        rank = matrix_rank(R)
        if rank == R.shape[0]:
            if rank == m - 1:
                break
            continue

        selection[n] = False

    # Find a linearly independent set of RR
    sigma = np.zeros((p - m + 1, p), dtype=int)

    for i, n in enumerate(reactions[~selection]):
        selection[n] = True

        _nu = nu[selection][:, inter][None]
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
    """Returns an array of possible response reaction routes for a given
    chemical formula array.

    Parameters
    ----------
    epsilon : ndarray (n, m)
        Chemical formula array of n elements by m molecular species.
    selection : ndarray (j,)
        Indices of the m species to be considered as terminal
    species : bool
        Return the indices of the chemical species used.

    Returns
    -------
    RER : ndarray (m, k)
        Possible response reaction routes.
    index : ndarray (j, k)
        Indices of the k terminal species use to produce the
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
