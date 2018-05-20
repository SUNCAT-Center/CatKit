import numpy as np
import networkx as nx


def raw_properties(
        atoms=None,
        atoms_parameters=None,
        connectivity=None):
    """Return all atom properties without manipulation."""
    fingerprint = np.concatenate(atoms_parameters)

    return fingerprint


def periodic_convolution(
        atoms=None,
        atoms_parameters=None,
        connectivity=None,
        d=0,
        normalize=False):
    """Return the square of each property with each atom. For distance 1 the
    covolutions returns the multiple of each property for all neighboring atom
    pairs.
    """
    if d == 0:
        convolution = np.dot(atoms_parameters,
                             atoms_parameters.T).diagonal()
    elif d == 1:
        V = np.dot(atoms_parameters, connectivity)
        convolution = np.dot(V, atoms_parameters.T).diagonal()
    else:
        raise ValueError('Only periodic_convolution d=0/1 supported')

    if normalize:
        convolution = np.sqrt(convolution / len(atoms))

    return convolution


def bonding_convolution(
        atoms=None,
        atoms_parameters=None,
        connectivity=None):
    """Perform convolution of metal atoms with bonded adsorbates."""
    # This is a CatKit convention
    bond_index = np.where(atoms.get_tags() == -1)[0]

    V = np.dot(atoms_parameters, connectivity)[:, bond_index]
    P = np.dot(V, atoms_parameters[:, bond_index].T).diagonal()
    convolution = P / connectivity[bond_index].sum()

    return convolution


def autocorrelation(
        atoms=None,
        atoms_parameters=None,
        connectivity=None,
        d=0):
    """Autocorrelation convolution for systems without pbc."""

    G = nx.Graph(connectivity)
    D = nx.floyd_warshall_numpy(G)
    S = np.zeros_like(D)
    S[D == d] = 1

    AC = np.dot(np.dot(atoms_parameters, S), atoms_parameters.T).diagonal()

    return AC
