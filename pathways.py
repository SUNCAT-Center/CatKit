import sqlite3
import numpy as np
import networkx as nx
from ase.data import atomic_numbers as an
import networkx.algorithms.isomorphism as iso
from ase.data import chemical_symbols as cs
import matplotlib.pyplot as plt
from rdkit.Chem.Draw import MolToFile
from rdkit.Chem import AllChem as Chem


def isomorphic_molecules(molecule1, molecule2):
    em = iso.numerical_edge_match('bonds', 1)
    nm = iso.numerical_node_match('atomic_number', 1)
    isomorphic = nx.is_isomorphic(
        molecule1,
        molecule2,
        edge_match=em,
        node_match=nm)

    return isomorphic


def get_bonds(G, u):

    bonds = 0
    for v in G.neighbors(u):
        bonds += G.get_edge_data(u, v)['bonds']

    return bonds


def plot_molecule(molecule, file_name=None):

    rdkG = get_rdkit_graph(molecule)
    MolToFile(rdkG, file_name)


def get_chemical_descriptor(graph, rank=1):

    symbols = nx.get_node_attributes(graph, 'symbol')
    chemical_tag = ''.join(np.sort(list(symbols.values())))

    if rank == 1:

        bonds = ''
        for u, v, data in graph.edges_iter(data=True):
            bond = symbols[u] + symbols[v]
            bonds += bond * data.get('bonds')

        bonds_tag = ''.join(sorted(bonds))

        return chemical_tag, bonds_tag

    return chemical_tag


def get_rdkit_graph(molecule):

    rdkitmol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())

    for j, data in molecule.nodes_iter(data=True):
        rdAtom = Chem.rdchem.Atom(data['symbol'])
        rdAtom.SetNumRadicalElectrons(int(data['radicals']))
        rdkitmol.AddAtom(rdAtom)

    rdBonds = Chem.rdchem.BondType
    orders = {
        '1': rdBonds.SINGLE,
        '2': rdBonds.DOUBLE,
        '3': rdBonds.TRIPLE}

    for u, v, data in molecule.edges_iter(data=True):

        order = orders[str(data['bonds'])]
        rdkitmol.AddBond(int(u), int(v), order)

    rdkitmol = rdkitmol.GetMol()
    rdkitmol = Chem.RemoveHs(rdkitmol)

    return rdkitmol


def get_smiles(molecule):

    rdkG = get_rdkit_graph(molecule)

    return Chem.MolToSmiles(rdkG)


class ReactionNetwork():
    """ A class for accessing a temporary SQLite database. This
    function works as a context manager and should be used as follows:

    with ReactionNetwork() as molnet:
        (Perform operation here)

    This syntax will automatically construct the temporary database,
    or access an existing one. Upon exiting the indentation, the
    changes to the database will be automatically committed.
    """

    def __init__(
            self,
            db_name='reaction-network.db'):
        """ The __init__ function is automatically called when the
        class is referenced.

        Args:
            db_name (str): Name of the database file to access. Will
            connect to 'molecule-network.db' by default.
        """

        self.db_name = db_name
        self.base_free_radicals = {
            'H': 1,
            'C': 4,
            'O': 2,
            'N': 3,
            'S': 2}

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
        product1 INT NOT NULL,
        product2 INT,
        FOREIGN KEY(reactant1) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(reactant2) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(product1) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(product2) REFERENCES molecules(molecule_pid),
        UNIQUE (reactant1, reactant2, product1, product2)
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS atoms(
        atom_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INT NOT NULL,
        node_id INT NOT NULL,
        atom_num INT NOT NULL,
        symbol CHAR NOT NULL,
        radicals INT NOT NULL,
        FOREIGN KEY(molecule_id) REFERENCES molecules(molecule_pid)
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS bonds(
        bond_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INT NOT NULL,
        node_id1 INT NOT NULL,
        node_id2 INT NOT NULL,
        nbonds INT NOT NULL,
        FOREIGN KEY(molecule_id) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(node_id1) REFERENCES atoms(atom_pid),
        FOREIGN KEY(node_id2) REFERENCES atoms(atom_pid)
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS positions(
        atom_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INT NOT NULL,
        atom_id INT NOT NULL,
        x_coord REAL NOT NULL,
        y_coord REAL NOT NULL,
        z_coord REAL NOT NULL,
        radicals INT NOT NULL,
        FOREIGN KEY(molecule_id) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(node_id1) REFERENCES atoms(atom_pid)
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS energies(
        atom_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INT NOT NULL,
        energy REAL NOT NULL,
        calculator CHAR,
        FOREIGN KEY(molecule_id) REFERENCES molecules(molecule_pid)
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
            molecule = nx.Graph()
            molecule.add_node(
                0,
                symbol=el,
                atomic_number=an[el],
                radicals=self.base_free_radicals[el])

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
            radicals = nodes[base_node]['radicals']
            base_el = molecule.node[base_node]['symbol']

            if radicals <= 0:
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

                G.node[base_node]['radicals'] -= 1
                G.add_node(
                    node_index,
                    symbol=el,
                    atomic_number=an[el],
                    radicals=self.base_free_radicals[el] - 1)
                G.add_edge(base_node, node_index, bonds=1)

                chemical_tag = ''.join(sorted(symbols + el))
                bonds_tag = ''.join(sorted(bonds + base_el + el))

                isomorph_found = False
                if chemical_tag not in molecules:
                    molecules[chemical_tag] = {}

                if bonds_tag not in molecules[chemical_tag]:
                    molecules[chemical_tag][bonds_tag] = []

                else:
                    for G1 in molecules[chemical_tag][bonds_tag]:
                        if isomorphic_molecules(G, G1):
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

                    radicals_new = nodes[existing_node]['radicals']
                    el = molecule.node[existing_node]['symbol']

                    if radicals_new <= 0:
                        continue

                    # Oxygen shouldn't bond to itself.
                    if base_el == 'O' and el == 'O':
                        continue

                    G = molecule.copy()

                    G.node[base_node]['radicals'] -= 1
                    G.node[existing_node]['radicals'] -= 1
                    G[base_node][existing_node]['bonds'] += 1

                    chemical_tag = ''.join(sorted(symbols))
                    bonds_tag = ''.join(sorted(bonds + base_el + el))

                    isomorph_found = False
                    if chemical_tag not in molecules:
                        molecules[chemical_tag] = {}

                    if bonds_tag not in molecules[chemical_tag]:
                        molecules[chemical_tag][bonds_tag] = []

                    else:
                        for G1 in molecules[chemical_tag][bonds_tag]:
                            if isomorphic_molecules(G, G1):
                                isomorph_found = True
                                break

                    if not isomorph_found:
                        new_molecules += [G]
                        molecules[chemical_tag][bonds_tag] += [G]

        return new_molecules, molecules

    def path_search(
            self,
            reconfiguration=True):

        self.reconfiguration = reconfiguration

        molecules = self.load_molecules(binned=True)

        pathways = []
        for data in molecules.values():
            for molecule_list in data.values():
                for molecule in molecule_list:

                    pathways += self.get_addition_paths(molecule, molecules)

        self.save_pathways(pathways)

        if reconfiguration:
            pathways = np.array(pathways)
            re_pathways = self.get_reconfiguration_paths(molecules, pathways)

            self.save_pathways(re_pathways)

        return pathways

    def get_addition_paths(
            self,
            molecule,
            molecules):

        disjoints, pathways = [], []
        for u, v, data in molecule.edges_iter(data=True):

            cut_molecule = molecule.copy()

            cut_molecule[u][v]['bonds'] -= 1
            if not cut_molecule[u][v]['bonds']:
                cut_molecule.remove_edge(u, v)

            isomorph_found = False
            for disjoint in disjoints:
                if isomorphic_molecules(cut_molecule, disjoint):
                    isomorph_found = True
                    break

            if not isomorph_found:
                disjoints += [cut_molecule]

                product_index = cut_molecule.graph['index']
                pieces = list(nx.connected_component_subgraphs(cut_molecule))

                addition_pathway = np.array([
                    None,
                    None,
                    product_index,
                    None])
                for i, piece in enumerate(pieces):

                    chemical_tag, bonds_tag = get_chemical_descriptor(piece)

                    for reactant in molecules[chemical_tag][bonds_tag]:
                        if isomorphic_molecules(piece, reactant):
                            pindex = reactant.graph['index']
                            addition_pathway[i] = pindex
                            break

                pathways += [addition_pathway]

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

                    for u, v, data in molecule.edges_iter(data=True):
                        bonds = data['bonds']
                        self.c.execute("""INSERT INTO bonds
                        (molecule_id, node_id1, node_id2, nbonds)
                        VALUES(?, ?, ?, ?)""", (molecule_pid, u, v, bonds))

                    for node, data in molecule.nodes_iter(data=True):
                        number = data['atomic_number']
                        radicals = data['radicals']
                        self.c.execute("""INSERT INTO atoms
                        (molecule_id, node_id, atom_num, symbol, radicals)
                        VALUES(?, ?, ?, ?, ?)""",
                                       (molecule_pid,
                                        node,
                                        number,
                                        cs[number],
                                        radicals))

    def save_pathways(self, pathways):

        for R1, R2, P1, P2 in pathways:

            self.c.execute("""INSERT INTO reactions
            (product1, product2, reactant1, reactant2)
            VALUES(?, ?, ?, ?)""", (P1, P2, R1, R2))

    def load_molecules(self, binned=False):

        cmd = """SELECT m.molecule_pid,
                        m.chemical_tag,
                        m.bonds_tag,
                        nodes,
                        bonds
        FROM molecules m LEFT JOIN
        (
         SELECT molecule_id,
         GROUP_CONCAT(node_id || ',' || atom_num || ',' || radicals, ';')
         as nodes
         FROM atoms
         GROUP BY molecule_id
        ) a
          ON a.molecule_id = m.molecule_pid LEFT JOIN
        (
         SELECT molecule_id,
         GROUP_CONCAT(node_id1 || ',' || node_id2 || ',' || nbonds, ';')
         as bonds
         FROM bonds
         GROUP BY molecule_id
        ) b
          ON b.molecule_id = m.molecule_pid
        """

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        molecules = {}

        for index, chemical_tag, bonds_tag, node_data, edge_data in fetch:

            molecule = nx.Graph(index=index)

            node_data = np.array([_.split(',')
                                  for _ in node_data.split(';')],
                                 dtype=int)
            nodes = [(node, {
                'atomic_number': n,
                'symbol': cs[n],
                'radicals': radicals})
                     for node, n, radicals in node_data]
            molecule.add_nodes_from(nodes)

            if edge_data:
                edges = np.array([_.split(',')
                                  for _ in edge_data.split(';')],
                                 dtype=int)
                molecule.add_weighted_edges_from(edges, weight='bonds')

            if binned:
                if chemical_tag not in molecules:
                    molecules[chemical_tag] = {}

                if bonds_tag not in molecules[chemical_tag]:
                    molecules[chemical_tag][bonds_tag] = []

                molecules[chemical_tag][bonds_tag] += [molecule]
            else:
                molecules[index] = molecule

        return molecules

    def load_pathways(self):

        cmd = """SELECT
        reactant1,
        reactant2,
        product1,
        product2
        FROM reactions"""

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        pathways = []
        for path in fetch:
            pathways += [[int(_) for _ in path if _]]

        pathways = np.array(pathways)

        return pathways

    def plot_reaction_network(self, file_name=None):

        cmd = """SELECT molecule_pid, chemical_tag
        FROM molecules"""

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        molecules = []
        for molecule in fetch:
            molecules += [(molecule[0],
                           {'tag': molecule[1], 'id': molecule[0]})]

        pathways = []
        for path in self.load_pathways():
            for R in path[:2]:
                for P in path[2:]:
                    if not R and not P:
                        pathways += [[R, P]]

        network = nx.Graph()
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

    def get_unsaturated_nodes(self, molecule, screen=None):

        unsaturated = []
        for node, data in molecule.nodes_iter(data=True):
            radicals = data['radicals']

            if screen in data:
                continue

            if radicals > 0:
                unsaturated += [node]

        return np.array(unsaturated)

    def get_cyclical_groups(self, molecule):
        """ Return a list of nodes that make up cyclical
        groups present in a molecule.
        """

        return nx.cycle_basis(molecule)

    def get_chirality(self, molecule):

        symbols = nx.get_node_attributes(molecule, 'symbol')

        for node, symbol in symbols.items():
            if symbol != 'C':
                continue

            neighbor_count = len(molecule.neighbors(node))

            if neighbor_count != 4:
                continue

            # Chirality is difficult to define correctly for
            # all molecules which contain multiple cycles.
            cycles = len(self.get_cyclical_groups(molecule))
            if cycles > 1:
                continue

            broken_molecule = molecule.copy()
            broken_molecule.remove_node(node)
            pieces = list(nx.connected_component_subgraphs(broken_molecule))
            n_pieces = len(pieces)

            chiral = True
            for i in range(n_pieces):
                piece = pieces.pop()
                for comparison_piece in pieces:

                    if isomorphic_molecules(
                            piece,
                            comparison_piece):
                        chiral = False
                        break

            if chiral:
                molecule.node[node]['chiral'] = 1

        return molecule

    def get_isomer_bonds(self, molecule):
        # Function is incomplete

        symbols = nx.get_node_attributes(molecule, 'symbol')

        for u, v, data in molecule.edges_iter(data=True):
            if data['bonds'] == 2:
                if symbols[u] == 'C' and symbols[v] == 'C':

                    cycles = len(self.get_cyclical_groups(molecule))
                    if cycles > 0:
                        continue

                    vbroken_molecule = molecule.copy()
                    vbroken_molecule.remove_node(u)

                    for node in vbroken_molecule.node_iter():
                        if not False:  # Check if node is connected to v
                            vbroken_molecule.remove_node()
                    vbroken_molecule.remove_node(v)

                    ubroken_molecule = molecule.copy()
                    ubroken_molecule.remove_node(v)

                    for node in ubroken_molecule.node_iter():
                        if not False:  # Check if node is connected to u
                            vbroken_molecule.remove_node(u)
                    vbroken_molecule.remove_node(u)

                    vpieces = list(nx.connected_component_subgraphs(
                        vbroken_molecule))
                    upieces = list(nx.connected_component_subgraphs(
                        ubroken_molecule))

                    if len(vpieces) == 2 and len(upieces) == 2:

                        isomer = False
                        if nx.is_isomorphic(
                                vpieces[0],
                                vpieces[1],
                                node_match=iso.numerical_node_match(
                                    'atomic_number', 1)):
                            isomer = False
                            break

                        if nx.is_isomorphic(
                                upieces[0],
                                upieces[1],
                                node_match=iso.numerical_node_match(
                                    'atomic_number', 1)):
                            isomer = False

                    if isomer:
                        # Store a bond as being chiral
                        molecule.node[node]['chiral'] = 1

    def get_reconfiguration_paths(
            self,
            molecules,
            pathways):

        ind_mol = self.load_molecules()

        new_pathways, reconfig = [], []
        for R1, R2, P1, P2 in pathways:
            if not R2 or P2:
                continue

            mol_R1 = ind_mol[R1]
            mol_R2 = ind_mol[R2]

            reconfigurations = [ind_mol[P1]]

            # Find potential bonding sites
            p1_bonds = self.get_unsaturated_nodes(
                ind_mol[R1])

            p2_bonds = self.get_unsaturated_nodes(
                ind_mol[R2]) + nx.number_of_nodes(ind_mol[R1])

            for b1 in p1_bonds:
                for b2 in p2_bonds:

                    Pt = nx.disjoint_union(mol_R1, mol_R2)
                    Pt.add_edge(b1, b2, bonds=1)

                    isomorph_found = False
                    for P in reconfigurations:
                        if isomorphic_molecules(Pt, P):
                            isomorph_found = True
                            break

                    if not isomorph_found:
                        reconfigurations += [Pt]

            del reconfigurations[0]
            reconfig_pathway = np.array([P1, None, None, None])
            for Pt in reconfigurations:
                chemical_tag, bonds_tag = get_chemical_descriptor(Pt)

                for P in molecules[chemical_tag][bonds_tag]:
                    if isomorphic_molecules(P, Pt):
                        reconfig_pathway[2] = P.graph['index']
                        pc = set([P.graph['index'], P1])
                        break

                if pc not in reconfig:
                    reconfig += [pc]
                    new_pathways += [reconfig_pathway]

        return new_pathways

#### NOT YET IMPLEMENTED ###
    
    # def get_substitution_paths(
    #         self,
    #         molecules,
    #         pathways):
    #     # Follows the form:
    #     # R1(-P2) + R2(-P2) --> P1 + P2

    #     ind_mol = self.load_molecules()

    #     for iP2, P2 in ind_mol.items():

    #         pbonds = self.get_unsaturated_nodes(P2)

    #         if not pbonds:
    #             continue

    #         sub_pathway = np.array([
                
                
                
    #         ])

    #         # Need to eliminate isomorphs to
    #         # not count double bonds
    #         joint_molecules = []
    #         for iR1, iR2, iP1, _ in pathways:

    #             R1 = ind_mol[iR1]

    #             for b in pbonds:
    #                 b2 = list(nx.get_node_attributes(
    #                     R,
    #                     'cut_node').keys())[0]
    #                 R_P2 = nx.disjoint_union(P2, R1)
    #                 R_P2.add_edge(b, b2, bonds=1)

    #                 isomorph_found = False
    #                 for G in joint_molecules:
    #                     if isomorphic_molecules(R_P2, G):
    #                         isomorph_found = True

    #                 if not isomorph_found:
    #                     joint_molecules += [R_P2]

    #             R2 = ind_mol[iR2]

    #             for b in pbonds:
    #                 b2 = list(nx.get_node_attributes(
    #                     R,
    #                     'cut_node').keys())[0]
    #                 R_P2 = nx.disjoint_union(P2, R2)
    #                 R_P2.add_edge(b, b2, bonds=1)

    #                 isomorph_found = False
    #                 for G in joint_molecules:
    #                     if isomorphic_molecules(R_P2, G):
    #                         isomorph_found = True

    #                 if not isomorph_found:
    #                     joint_molecules += [R_P2]


    #             for R in joint_molecules:
    #                 (chemical_tag,
    #                  bonds_tag) = get_chemical_descriptor(R)

    #                 for G in molecules[chemical_tag][bonds_tag]:
    #                     if isomorphic_molecules(R, G):

    #                         break

    #                 pathways += [sub_pathway]
