import sqlite3
from sqlite3 import IntegrityError
import numpy as np
import networkx as nx
from ase.data import atomic_numbers as an
import networkx.algorithms.isomorphism as iso
from ase.data import chemical_symbols as cs
import matplotlib.pyplot as plt


def get_chemical_tag(graph):

    attr = nx.get_node_attributes(graph, 'symbol')
    tag = ''.join(np.sort(list(attr.values())))

    return tag


class ReactionNetwork():
    """ A class for accessing a temporary SQLite database. This
    function works as a context manager and should be used as follows:

    with ReactionNetwork() as molnet:
        (Perform operation here)

    This syntax will automatically construct the temporary database,
    or access an existing one. Upon exiting the indentation, the
    changes to the database will be automatically commited.
    """

    def __init__(self, db_name='reaction-network.db'):
        """ The __init__ function is automatically called when the
        class is referenced.

        Args:
            db_name (str): Name of the database file to access. Will
            connect to 'molecule-network.db' by default.
        """

        self.db_name = db_name

    def __enter__(self):
        """ This function is automatically called whenever the class
        is used together with a 'with' statement.
        """

        self.con = sqlite3.connect(self.db_name)
        self.c = self.con.cursor()
        self.create_table()

        return self

    def __exit__(self, type, value, tb):
        """ Upon exiting the 'with' statement, __exit__ is called.
        """

        self.con.commit()
        self.con.close()

    def create_table(self):
        """ Creates the database table framework used in SQLite.
        This includes 3 tables: molecules, reasctions, and bonds.
        """

        self.c.execute("""CREATE TABLE IF NOT EXISTS molecules(
        molecule_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        chemical_tag CHAR NOT NULL
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS reactions(
        reaction_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        reactant1 INT NOT NULL,
        reactant2 INT NOT NULL,
        product INT NOT NULL,
        FOREIGN KEY(reactant1) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(reactant2) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(product) REFERENCES molecules(molecule_pid)
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS atoms(
        atom_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INT NOT NULL,
        node_id INT NOT NULL,
        atom_num INT NOT NULL,
        symbol CHAR NOT NULL,
        FOREIGN KEY(molecule_id) REFERENCES molecules(molecule_pid)
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS bonds(
        bond_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INT NOT NULL,
        node_id1 INT NOT NULL,
        node_id2 INT NOT NULL,
        FOREIGN KEY(molecule_id) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(node_id1) REFERENCES atoms(atom_pid),
        FOREIGN KEY(node_id2) REFERENCES atoms(atom_pid)
        )""")

    def molecule_search(
            self,
            element_pool={'C': 2, 'H': 6},
            max_bonds={'H': 1, 'C': 4, 'O': 2}):

        self.element_pool = element_pool
        self.max_bonds = max_bonds

        molecules = {}
        search_molecules = []
        for el in list(self.element_pool.keys()):
            molecule = nx.Graph()
            molecule.add_node(
                0,
                symbol=el,
                atomic_number=an[el],
                max_bonds=self.max_bonds[el])

            search_molecules += [molecule]
            molecules[el] = [molecule]

        for molecule in search_molecules:

            # Recusive use of molecules
            new_molecules, molecules = self.branch_molecule(
                molecule,
                molecules)

            search_molecules += new_molecules

        return molecules

    def branch_molecule(self, molecule, molecules):

        # In the future, making holding molecules in memory optional

        new_molecules = []
        nodes = molecule.node

        symbol_tags = list(nx.get_node_attributes(molecule, 'symbol').values())
        symbol_count = {sym: symbol_tags.count(
            sym) for sym in set(symbol_tags)}

        for base_node in nodes:

            # Check if an additional bond can be formed
            max_bonds = nodes[base_node]['max_bonds']
            edge_count = len(molecule.edges([base_node]))

            if edge_count >= max_bonds:
                continue

            for el, cnt in self.element_pool.items():

                # Skip elements which has been fully saturated.
                if symbol_count.get(el, 0) >= cnt:
                    continue

                # Oxygen shouldn't bond to itself.
                if nodes[base_node]['symbol'] == 'O' and el == 'O':
                    continue

                G = molecule.copy()

                node_index = G.number_of_nodes()

                G.add_node(
                    node_index,
                    symbol=el,
                    atomic_number=an[el],
                    max_bonds=self.max_bonds[el])

                G.add_edge(base_node, node_index)

                tag = ''.join(np.sort(symbol_tags + [el]))

                isomorph_found = False
                if tag not in list(molecules.keys()):
                    molecules[tag] = []

                else:
                    for G1 in molecules[tag]:
                        if nx.is_isomorphic(
                                G,
                                G1,
                                node_match=iso.numerical_node_match(
                                    'atomic_number', 1)):
                            isomorph_found = True
                            break

                if not isomorph_found:
                    new_molecules += [G]
                    molecules[tag] += [G]

        return new_molecules, molecules

    def path_search(self):

        molecules = self.load_molecules()

        pathways = []
        for tag, graphs in molecules.items():
            for molecule in graphs:

                pathways += self.branch_paths(molecule, molecules)

        return pathways

    def branch_paths(self, molecule, molecules):

        disjoints, pathways = [], []
        for edge in molecule.edges_iter():
            cut_molecule = molecule.copy()
            cut_molecule.remove_edge(*edge)

            isomorph_found = False
            for disjoint in disjoints:
                if nx.is_isomorphic(
                    cut_molecule,
                    disjoint,
                    node_match=iso.numerical_node_match(
                        'atomic_number',
                        1)):
                    isomorph_found = True
                    break

            if not isomorph_found:
                disjoints += [cut_molecule]

                product_index = cut_molecule.graph['index']
                pieces = nx.connected_component_subgraphs(cut_molecule)

                pathway = [product_index]
                for i, piece in enumerate(pieces):

                    tag = get_chemical_tag(piece)

                    for reactant in molecules[tag]:
                        if nx.is_isomorphic(
                            piece,
                            reactant,
                            node_match=iso.numerical_node_match(
                                'atomic_number',
                                1)):

                            reactant_index = reactant.graph['index']
                            pathway += [reactant_index]
                            break

                pathways += [pathway]

        return pathways

    def save_molecules(self, molecules):

        for tag, molecule_list in molecules.items():
            for molecule in molecule_list:

                self.c.execute("""INSERT INTO molecules (chemical_tag)
                VALUES(?)""", (tag,))

                self.c.execute("""SELECT last_insert_rowid()""")
                molecule_pid = self.c.fetchone()[0]

                for edge in molecule.edges_iter():
                    n1 = edge[0]
                    n2 = edge[1]
                    self.c.execute("""INSERT INTO bonds
                    (molecule_id, node_id1, node_id2)
                    VALUES(?, ?, ?)""", (molecule_pid, n1, n2))

                numbers = nx.get_node_attributes(molecule, 'atomic_number')
                for node, number in numbers.items():
                    self.c.execute("""INSERT INTO atoms
                    (molecule_id, node_id, atom_num, symbol)
                    VALUES(?, ?, ?, ?)""", (molecule_pid, node, number, cs[number]))

    def save_pathways(self, pathways):

        for product, reactant1, reactant2 in pathways:

            self.c.execute("""INSERT INTO reactions
            (product, reactant1, reactant2)
            VALUES(?, ?, ?)""", (product, reactant1, reactant2))

    def load_molecules(self):

        cmd = """SELECT molecule_id, chemical_tag,
        GROUP_CONCAT(nodes, ';'), GROUP_CONCAT(bonds, ';')
        FROM (SELECT a.molecule_id,  m.chemical_tag,
        a.node_id || ',' || a.atom_num AS nodes,
        b.node_id1 || ',' || b.node_id2 AS bonds
        FROM atoms a
        JOIN molecules m on a.molecule_id = m.molecule_pid
        LEFT OUTER JOIN bonds b on a.molecule_id = b.molecule_id)
        GROUP BY molecule_id"""

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        molecules = {}
        for index, tag, node_data, edge_data in fetch:

            molecule = nx.Graph(index=index)

            node_data = np.array([_.split(',')
                                  for _ in node_data.split(';')], dtype=int)
            nodes = [(node, {'atomic_number': n, 'symbol': cs[n]})
                     for node, n in node_data]
            molecule.add_nodes_from(nodes)

            if edge_data:
                edges = np.array([_.split(',')
                                  for _ in edge_data.split(';')], dtype=int)
                molecule.add_edges_from(edges)

            if tag not in list(molecules.keys()):
                molecules[tag] = []

            molecules[tag] += [molecule]

        return molecules

    def plot_molecule(self, molecule, file_name=None):

        colors = {'H': 'w', 'C': 'grey', 'O': 'r'}

        symbols = nx.get_node_attributes(molecule, 'symbol')
        color = [colors[sym] for _, sym in symbols.items()]

        plt.figure(figsize=(3, 3))
        nx.draw_networkx(
            molecule,
            labels=nx.get_node_attributes(molecule, 'symbol'),
            node_color=color)
        plt.draw()
        plt.axis('off')

        if file_name:
            plt.savefig(file_name)

    def plot_reaction_network(self, file_name=None):

        cmd = """SELECT molecule_pid, chemical_tag
        FROM molecules"""

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        molecules = []
        for molecule in fetch:
            molecules += [(molecule[0], {'tag': molecule[1]})]

        cmd = """SELECT
        product || ',' || reactant1  AS path1,
        product || ',' || reactant2  AS path2
        FROM reactions"""

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        pathways = []
        for paths in fetch:
            pathways += [path.split(',') for path in paths]

        pathways = np.array(pathways, dtype=int)

        network = nx.Graph()
        network.add_nodes_from(molecules)
        network.add_edges_from(pathways)

        plt.figure(figsize=(6, 6))
        nx.draw_networkx(
            network,
            labels=nx.get_node_attributes(network, 'tag'))
        plt.draw()
        plt.axis('off')

        if file_name:
            plt.savefig(file_name)
