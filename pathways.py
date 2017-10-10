import sqlite3
import numpy as np
import networkx as nx
from ase.data import atomic_numbers as an
import networkx.algorithms.isomorphism as iso
from ase.data import chemical_symbols as cs
import matplotlib.pyplot as plt


def plot_molecule(molecule, file_name=None):

    colors = {'H': 'w', 'C': 'grey', 'O': 'r'}

    symbols = nx.get_node_attributes(molecule, 'symbol')
    color = [colors[sym] for _, sym in list(symbols.items())]

    fig, ax = plt.subplots(figsize=(2, 2))
    nx.draw_networkx(
        molecule,
        labels=nx.get_node_attributes(molecule, 'symbol'),
        node_color=color)
    plt.draw()

    plt.title(molecule.graph['index'])
    plt.axis('off')

    if file_name:
        plt.savefig(file_name)

    plt.close()


def get_chemical_descriptor(graph, rank=1):

    symbols = nx.get_node_attributes(graph, 'symbol')
    chemical_tag = ''.join(np.sort(list(symbols.values())))

    if rank == 1:
        edge_list = list(graph.edges())

        if edge_list:
            bonds = ''
            for e1, e2 in edge_list:
                bonds += symbols[e1] + symbols[e2]
            bonds_tag = ''.join(sorted(bonds))
        else:
            bonds_tag = ''

        return chemical_tag, bonds_tag

    return chemical_tag


class ReactionNetwork():
    """ A class for accessing a temporary SQLite database. This
    function works as a context manager and should be used as follows:

    with ReactionNetwork() as molnet:
        (Perform operation here)

    This syntax will automatically construct the temporary database,
    or access an existing one. Upon exiting the indentation, the
    changes to the database will be automatically commited.
    """

    def __init__(
            self,
            db_name='reaction-network.db',
            max_bonds={'H': 1, 'C': 4, 'O': 2}):
        """ The __init__ function is automatically called when the
        class is referenced.

        Args:
            db_name (str): Name of the database file to access. Will
            connect to 'molecule-network.db' by default.
        """

        self.db_name = db_name
        self.max_bonds = max_bonds

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
        This includes 4 tables: molecules, reasctions, atoms, and bonds.
        """

        self.c.execute("""CREATE TABLE IF NOT EXISTS molecules(
        molecule_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        chemical_tag CHAR NOT NULL,
        bonds_tag CHAR NOT NULL
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS reactions(
        reaction_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        reactant1 INT NOT NULL,
        reactant2 INT,
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
            multiple_bond_search=False):

        self.multiple_bond_search = multiple_bond_search
        self.element_pool = element_pool

        molecules = {}
        search_molecules = []
        for el in self.element_pool:
            molecule = nx.MultiGraph()
            molecule.add_node(
                0,
                symbol=el,
                atomic_number=an[el],
                max_bonds=self.max_bonds[el])

            search_molecules += [molecule]
            molecules[el] = {'': [molecule]}

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

        symbols, bonds = get_chemical_descriptor(molecule)
        symbol_count = {sym: symbols.count(sym)
                        for sym in set(symbols)}

        for base_node in nodes:

            # Check if an additional bond can be formed
            max_bonds = nodes[base_node]['max_bonds']
            edge_count = len(molecule.edges([base_node]))
            base_el = molecule.node[base_node]['symbol']

            if edge_count >= max_bonds:
                continue

            # Creating new nodes
            for el, cnt in list(self.element_pool.items()):

                # Skip elements exceed specified limit
                if symbol_count.get(el, 0) >= cnt:
                    continue

                # Oxygen shouldn't bond to itself.
                if base_el == 'O' and el == 'O':
                    continue

                G = molecule.copy()

                node_index = G.number_of_nodes()

                G.add_node(
                    node_index,
                    symbol=el,
                    atomic_number=an[el],
                    max_bonds=self.max_bonds[el])

                G.add_edge(base_node, node_index)

                chemical_tag = ''.join(sorted(symbols + el))
                bonds_tag = ''.join(sorted(bonds + base_el + el))

                isomorph_found = False
                if chemical_tag not in molecules:
                    molecules[chemical_tag] = {}

                if bonds_tag not in molecules[chemical_tag]:
                    molecules[chemical_tag][bonds_tag] = []

                else:
                    for G1 in molecules[chemical_tag][bonds_tag]:
                        if nx.is_isomorphic(
                                G,
                                G1,
                                node_match=iso.numerical_node_match(
                                    'atomic_number', 1)):
                            isomorph_found = True
                            break

                if not isomorph_found:
                    new_molecules += [G]
                    molecules[chemical_tag][bonds_tag] += [G]

            # Bonds to existing nodes
            if self.multiple_bond_search:
                for existing_node in nodes:

                    # Atoms don't bond to themselves
                    if existing_node == base_node:
                        continue

                    max_bonds_new = nodes[existing_node]['max_bonds']
                    edge_count_new = len(molecule.edges([existing_node]))
                    el = molecule.node[existing_node]['symbol']

                    if edge_count_new >= max_bonds_new:
                        continue

                    # Oxygen shouldn't bond to itself.
                    if base_el == 'O' and el == 'O':
                        continue

                    G = molecule.copy()

                    G.add_edge(base_node, existing_node)

                    chemical_tag = ''.join(sorted(symbols))
                    bonds_tag = ''.join(sorted(bonds + base_el + el))

                    isomorph_found = False
                    if chemical_tag not in molecules:
                        molecules[chemical_tag] = {}

                    if bonds_tag not in molecules[chemical_tag]:
                        molecules[chemical_tag][bonds_tag] = []

                    else:
                        for G1 in molecules[chemical_tag][bonds_tag]:
                            if nx.is_isomorphic(
                                    G,
                                    G1,
                                    node_match=iso.numerical_node_match(
                                        'atomic_number', 1)):
                                isomorph_found = True
                                break

                    if not isomorph_found:
                        new_molecules += [G]
                        molecules[chemical_tag][bonds_tag] += [G]

        return new_molecules, molecules

    def path_search(self):

        molecules = self.load_molecules()

        pathways = []
        for chemical_tag, data in molecules.items():
            for bonds_tag, molecule_list in data.items():
                for molecule in molecule_list:

                    pathways += self.branch_paths(molecule, molecules)

        return pathways

    def branch_paths(self, molecule, molecules, substitution=True):

        disjoints, pathways = [], []
        for u, v in molecule.edges():
            cut_molecule = molecule.copy()
            cut_molecule.remove_edge(u, v)

            isomorph_found = False
            for disjoint in disjoints:
                if nx.is_isomorphic(
                        cut_molecule,
                        disjoint,
                        node_match=iso.numerical_node_match(
                            'atomic_number', 1)):
                    isomorph_found = True
                    break

            if not isomorph_found:
                disjoints += [cut_molecule]

                product_index = cut_molecule.graph['index']
                pieces = nx.connected_component_subgraphs(cut_molecule)
                pieces_index = []

                pathway = [product_index]
                for piece in pieces:

                    chemical_tag, bonds_tag = get_chemical_descriptor(piece)

                    for reactant in molecules[chemical_tag][bonds_tag]:
                        if nx.is_isomorphic(
                                piece,
                                reactant,
                                node_match=iso.numerical_node_match(
                                    'atomic_number', 1)):
                            pindex = reactant.graph['index']
                            pieces_index += [pindex]
                            reactant_index = pindex
                            pathway += [reactant_index]
                            break

                pathways += [pathway]

        return pathways

    def save_molecules(self, molecules):

        for chemical_tag, data in list(molecules.items()):
            for bonds_tag, molecule_list in list(data.items()):
                for molecule in molecule_list:

                    self.c.execute("""INSERT INTO molecules
                    (chemical_tag, bonds_tag)
                    VALUES(?, ?)""", (chemical_tag, bonds_tag))

                    self.c.execute("""SELECT last_insert_rowid()""")
                    molecule_pid = self.c.fetchone()[0]

                    for edge in molecule.edges():
                        n1 = edge[0]
                        n2 = edge[1]
                        self.c.execute("""INSERT INTO bonds
                        (molecule_id, node_id1, node_id2)
                        VALUES(?, ?, ?)""", (molecule_pid, n1, n2))

                    numbers = nx.get_node_attributes(molecule, 'atomic_number')
                    for node, number in list(numbers.items()):
                        self.c.execute("""INSERT INTO atoms
                        (molecule_id, node_id, atom_num, symbol)
                        VALUES(?, ?, ?, ?)""",
                                       (molecule_pid,
                                        node,
                                        number,
                                        cs[number]))

    def save_pathways(self, pathways):

        for pathway in pathways:
            if len(pathway) == 2:
                product, reactant1 = pathway
                reactant2 = None
            elif len(pathway) == 3:
                product, reactant1, reactant2 = pathway

            self.c.execute("""INSERT INTO reactions
            (product, reactant1, reactant2)
            VALUES(?, ?, ?)""", (product, reactant1, reactant2))

    def load_molecules(self):

        cmd = """SELECT m.molecule_pid,
                        m.chemical_tag,
                        m.bonds_tag,
                        nodes,
                        bonds
        FROM molecules m LEFT JOIN
        (
         SELECT molecule_id,
         GROUP_CONCAT(node_id || ',' || atom_num, ';') as nodes
         FROM atoms
         GROUP BY molecule_id
        ) a
          ON a.molecule_id = m.molecule_pid LEFT JOIN
        (
         SELECT molecule_id,
         GROUP_CONCAT(node_id1 || ',' || node_id2, ';') as bonds
         FROM bonds
         GROUP BY molecule_id
        ) b
          ON b.molecule_id = m.molecule_pid
        """

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        molecules = {}
        for index, chemical_tag, bonds_tag, node_data, edge_data in fetch:

            molecule = nx.MultiGraph(index=index)

            node_data = np.array([_.split(',')
                                  for _ in node_data.split(';')],
                                 dtype=int)
            nodes = [(node, {'atomic_number': n, 'symbol': cs[n]})
                     for node, n in node_data]
            molecule.add_nodes_from(nodes)

            if edge_data:
                edges = np.array([_.split(',')
                                  for _ in edge_data.split(';')],
                                 dtype=int)
                molecule.add_edges_from(edges)

            if chemical_tag not in molecules:
                molecules[chemical_tag] = {}

            if bonds_tag not in molecules[chemical_tag]:
                molecules[chemical_tag][bonds_tag] = []

            molecules[chemical_tag][bonds_tag] += [molecule]

        return molecules

    def plot_reaction_network(self, file_name=None):

        cmd = """SELECT molecule_pid, chemical_tag
        FROM molecules"""

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        molecules = []
        for molecule in fetch:
            molecules += [(molecule[0],
                           {'tag': molecule[1], 'id': molecule[0]})]

        cmd = """SELECT
        product || ',' || reactant1  AS path1,
        product || ',' || reactant2  AS path2
        FROM reactions"""

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        pathways = []
        for paths in fetch:
            pathways += [path.split(',') for path in paths if path is not None]

        pathways = np.array(pathways, dtype=int)

        network = nx.MultiGraph()
        network.add_nodes_from(molecules)
        network.add_edges_from(pathways)

        plt.figure(figsize=(6, 6))
        nx.draw_networkx(
            network,
            labels=nx.get_node_attributes(network, 'id'))
        plt.draw()
        plt.axis('off')

        if file_name:
            plt.savefig(file_name)

        plt.close()
