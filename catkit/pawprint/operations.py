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


def layered_sum(
        atoms=None,
        atoms_parameters=None,
        connectivity=None):
    """Sum of the properties in a layer as indicated by catkit tags."""
    tags = atoms.get_tags()
    tags -= min(tags)
    LS = np.array([np.bincount(tags, weights=ap) for ap in atoms_parameters])
    LS = LS[LS != 0]

    return LS


def local_ads_metal_fp(
        atoms=None,
        atoms_parameters=None,
        connectivity=None,
        fuse=False):
    """Sum of the differences in properties of the atoms in the
       metal-adsorbate interface

       Parameters
       ----------
       atoms : ase Atoms or catkit gratoms object.
       atoms_parameters : ndarray(n, )
           a list of chemical properties to construct the fingerprints.
       connectivity : ndarray (n,)
           Connectivity of the adsorption sites
       weigthed : boolean
           fingerprints are weightd by the stoichiometric ratio of
           the bimetals.
    """
    bond_index = np.where(atoms.get_tags() == -1)[0]
    fp = np.empty([len(atoms_parameters), len(bond_index)])
    for i, bi in enumerate(bond_index):
        bonded_ap = atoms_parameters[:, np.where(connectivity[bi] == 1)[0]]
        fp[:, i] = np.mean(bonded_ap -
                atoms_parameters[:, bi].reshape(-1, 1), axis=1)
    if not fuse:
        return fp.reshape(-1)
    else:
       return fp.sum(axis=1)

def bimetal_fp(
        atoms=None,
        atoms_parameters=None,
        connectivity=None):
    """The differences in properties of the atoms in the
       metal-adsorbate interface

       Parameters
       ----------
       atoms : ase Atoms or catkit gratoms object.
       atoms_parameters : ndarray(n, )
           a list of chemical properties to construct the fingerprints.
    """
    fp = np.zeros([len(atoms_parameters)])
    metals = set(atoms.get_chemical_symbols())
    uap = []
    for ap in atoms_parameters:
        if not  np.isnan(ap).any():
            uap += [np.unique(np.round(ap, 3)).tolist()]
        else:
            uap += [[np.nan, np.nan]]
    for i, ap in enumerate(uap):
        if len(ap) == 1:
            uap[i] = [ap[0], ap[0]]

    if len(metals) == 1:
        fp = fp.reshape(-1)
    elif len(metals) == 2:
        fp = np.diff(uap).reshape(-1)
    else:
        raise NotImplementedError("""This operation is restricted to single
                                  and binary metal system only.""")
    return fp


def derived_fp(
        atoms=None,
        atoms_parameters=None,
        connectivity=None,
        fp_1=None,
        fp_2=None,
        n_1=None,
        n_2=None,
        op=None):
    """
    NOTE : This is a work in progress. I'll redesign the whole thing to allow
           for arithmetic manipulation of two fingerprints.
       Given two fingerprints vector, it will perform arithmetic operation to
       design new fingerprints.

       Available operations:
           add : adds two fingerprints of equal length raised to their given
                 power.
           subtract : subtracts two fingerprints of equal length raised to
                      their given power.
           mulltiply : multiply two fingerprints of equal length raised to
                       their given power.
           divide : divide two fingerprints of equal length raised to their
                    given power."""
    if op == 'add':
        return fp_1 ** n_1 + fp_2 ** n_2
    elif op == 'subtract':
        return fp_1 ** n_1 - fp_2 ** n_2
    elif op == 'divide':
        return fp_1 ** n_1 / fp_2 ** n_2
    elif op == 'multiply':
        return fp_1 ** n_1 * fp_2 ** n_2
