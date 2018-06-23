import networkx.algorithms.isomorphism as iso
import networkx as nx
import numpy as np


def connectivity_to_edges(connectivity, indices=None):
    """Convert a Numpy connectivity matrix into a list of NetworkX
    compatible edges.
    """
    if indices is None:
        indices = np.arange(connectivity.shape[0], dtype=int)

    edges = []
    for i, c in enumerate(connectivity):
        lower_diagonal = c[:i + 1]

        for j, v in enumerate(lower_diagonal):
            edges += [(indices[i], indices[j], 1)] * int(v)
    return edges


def isomorphic_molecules(graph0, graph1):
    """Check whether two molecule graphs are isomorphic."""
    em = iso.numerical_edge_match('bonds', 1)
    nm = iso.numerical_node_match('number', 1)

    isomorphic = nx.is_isomorphic(graph0, graph1, edge_match=em, node_match=nm)

    return isomorphic
