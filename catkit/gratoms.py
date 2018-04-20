from builtins import super
import networkx as nx
from networkx import Graph, MultiGraph
from ase import Atoms, Atom
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

    def __init__(self,
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
                 edges=None):
        super().__init__(symbols, positions, numbers, tags, momenta, masses,
                         magmoms, charges, scaled_positions, cell, pbc,
                         celldisp, constraint, calculator, info)

        if isinstance(edges, np.ndarray):
            if self.pbc.any():
                self._graph = MultiGraph(edges)
            else:
                self._graph = Graph(edges)
        else:
            if self.pbc.any():
                self._graph = MultiGraph()
            else:
                self._graph = Graph()

        nodes = [[i, {
            'number': n
        }] for i, n in enumerate(self.arrays['numbers'])]
        self._graph.add_nodes_from(nodes)

        if isinstance(edges, list):
            self._graph.add_edges_from(edges)

        self._surface_atoms = None

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

    @property
    def degree(self):
        degree = self._graph.degree
        return np.array([_[1] for _ in degree])

    @property
    def connectivity(self):
        connectivity = nx.to_numpy_matrix(self._graph).astype(int)
        return connectivity

    def get_surface_atoms(self):
        """Return surface atoms."""
        return self._surface_atoms

    def set_surface_atoms(self, surface_atoms):
        """Assign surface atoms."""
        self._surface_atoms = surface_atoms

    def get_neighbor_symbols(self, u):
        """Get chemical symbols for neighboring atoms of u."""
        neighbors = list(self._graph[u])

        return sym[self.arrays['numbers'][neighbors]]

    def is_isomorph(self, other):
        """Check if isomorphic by bond count and atomic number."""
        isomorphic = nx.is_isomorphic(
            self._graph, other._graph, edge_match=em, node_match=nm)

        return isomorphic

    def get_chemical_tags(self, rank=2):
        """Generate a hash descriptive of the chemical formula (rank 0)
        or include bonding (rank 1).
        """
        cnt = np.bincount(self.arrays['numbers'])
        composition = ','.join(cnt.astype(str))

        if rank == 1:
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
            a = np.zeros((n1 + n2, ) + a1.shape[1:], a1.dtype)
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
            a = np.empty((n1 + n2, ) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            if name == 'masses':
                a[:n1] = self.get_masses()[:n1]
            else:
                a[:n1] = 0

            self.set_array(name, a)

        if isinstance(other, Gratoms):
            if isinstance(self._graph, nx.MultiGraph) & \
               isinstance(other._graph, nx.Graph):
                other._graph = nx.MultiGraph(other._graph)

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

        n = len(self)
        i = np.arange(n)[i]

        if len(self._constraints) > 0:
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

    def __imul__(self, m):
        """In-place repeat of atoms."""
        if isinstance(m, int):
            m = (m, m, m)

        for x, vec in zip(m, self._cell):
            if x != 1 and not vec.any():
                raise ValueError(
                    'Cannot repeat along undefined lattice vector')

        if self.pbc.any() and len(self.edges()) > 0:
            raise ValueError("Edge conservation not currently supported with "
                             "pbc. Remove pbc or edges first.")

        M = np.product(m)
        n = len(self)

        for name, a in self.arrays.items():
            self.arrays[name] = np.tile(a, (M, ) + (1, ) * (len(a.shape) - 1))
            cgraph = self._graph.copy()

        positions = self.arrays['positions']
        i0 = 0

        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), self._cell)
                    i0 = i1
                    if m0 + m1 + m2 != 0:
                        self._graph = nx.disjoint_union(self._graph, cgraph)

        if self.constraints is not None:
            self.constraints = [c.repeat(m, n) for c in self.constraints]

        self._cell = np.array([m[c] * self._cell[c] for c in range(3)])

        return self
