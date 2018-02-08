from builtins import super
from ase import Atoms, Atom
import networkx as nx
from networkx import Graph, MultiGraph
from ase.data import chemical_symbols
import networkx.algorithms.isomorphism as iso
import numpy as np
import copy
sym = np.array(chemical_symbols)
em = iso.numerical_edge_match('bonds', 1)
nm = iso.numerical_node_match('number', 1)


class Gratoms(Atoms):
    """Graph based atoms object.

    An Integrated class for an ASE atoms object with a corresponding
    Networkx Graph.
    """

    def __init__(
            self,
            symbols=None,
            positions=None,
            numbers=None,
            tags=None,
            momenta=None,
            masses=None,
            magmoms=None,
            charges=None,
            scaled_positions=None,
            cell=None,
            pbc=None,
            celldisp=None,
            constraint=None,
            calculator=None,
            info=None,
            edges=None
    ):
        super().__init__(
            symbols, positions, numbers, tags, momenta, masses,
            magmoms, charges, scaled_positions, cell, pbc, celldisp,
            constraint, calculator, info)

        if self.pbc.any():
            self._graph = MultiGraph()
        else:
            self._graph = Graph()
        if edges:
            self._graph.add_edges_from(edges, bonds=1)
        nodes = [[i, {'number': n}]
                 for i, n in enumerate(self.arrays['numbers'])]
        self._graph.add_nodes_from(nodes)

    @property
    def graph(self):
        return self._graph

    @property
    def nodes(self):
        return self._graph.nodes

    @property
    def edges(self):
        return self._graph.edges

    @property
    def adj(self):
        return self._graph.adj

    def get_neighbor_symbols(self, u):
        """ Get chemical symbols for neighboring atoms of u."""
        neighbors = list(self._graph[u])

        return sym[self.arrays['numbers'][neighbors]]

    def is_isomorph(self, other):
        """Check if isomorphic by bond count and atomic number."""
        isomorphic = nx.is_isomorphic(
            self._graph,
            other._graph,
            edge_match=em,
            node_match=nm)

        return isomorphic

    def get_chemical_tags(self, rank=1):
        """Generate a hash descriptive of the chemical formula (rank 0)
        or include bonding (rank 1).
        """
        cnt = np.bincount(self.arrays['numbers'])
        composition = ','.join(cnt.astype(str))

        if rank != 1:
            return composition[2:]

        for adj in self.adj.items():

            num = self.arrays['numbers'][list(adj[1].keys())]
            cnt += np.bincount(num, minlength=len(cnt))

        bonding = ','.join(cnt.astype(str))

        return composition[2:], bonding[2:]

    def get_unsaturated_nodes(self, screen=None):

        unsaturated = []
        for node, data in self.nodes(data=True):
            radicals = data['valence']

            if screen in data:
                continue

            if radicals > 0:
                unsaturated += [node]

        return np.array(unsaturated)

    def copy(self):
        """Return a copy."""
        atoms = self.__class__(cell=self._cell, pbc=self._pbc, info=self.info)

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a.copy()
        atoms.constraints = copy.deepcopy(self.constraints)
        atoms._graph = self.graph.copy()

        return atoms

    def __iadd__(self, other):
        """Extend atoms object by appending atoms from *other*."""
        if isinstance(other, Atom):
            other = self.__class__([other])

        n1 = len(self)
        n2 = len(other)

        for name, a1 in self.arrays.items():
            a = np.zeros((n1 + n2,) + a1.shape[1:], a1.dtype)
            a[:n1] = a1
            if name == 'masses':
                a2 = other.get_masses()
            else:
                a2 = other.arrays.get(name)
            if a2 is not None:
                a[n1:] = a2
            self.arrays[name] = a

        for name, a2 in other.arrays.items():
            if name in self.arrays:
                continue
            a = np.empty((n1 + n2,) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            if name == 'masses':
                a[:n1] = self.get_masses()[:n1]
            else:
                a[:n1] = 0

            self.set_array(name, a)

        self._graph = nx.disjoint_union(self._graph, other._graph)

        return self

    def __delitem__(self, i):
        from ase.constraints import FixAtoms
        for c in self._constraints:
            if not isinstance(c, FixAtoms):
                raise RuntimeError('Remove constraint using set_constraint() '
                                   'before deleting atoms.')

        if isinstance(i, (list, int)):
            # Make sure a list of booleans will work correctly and not be
            # interpreted at 0 and 1 indices.
            i = np.atleast_1d(i)

        if len(self._constraints) > 0:
            n = len(self)
            i = np.arange(n)[i]
            if isinstance(i, int):
                i = [i]
            constraints = []
            for c in self._constraints:
                c = c.delete_atoms(i, n)
                if c is not None:
                    constraints.append(c)
            self.constraints = constraints

        mask = np.ones(len(self), bool)
        mask[i] = False

        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]

        self._graph.remove_nodes_from(i)
        mapping = dict(zip(np.where(mask)[0], np.arange(len(self))))
        nx.relabel_nodes(self._graph, mapping, copy=False)
