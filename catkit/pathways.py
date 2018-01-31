import sqlite3
import numpy as np
import networkx as nx
from ase.data import atomic_numbers as an
import networkx.algorithms.isomorphism as iso
from ase.data import chemical_symbols as cs
import matplotlib.pyplot as plt
from ase import Atoms
try:
    from rdkit.Chem.Draw import MolToFile
    from rdkit.Chem import AllChem as Chem
except(ImportError):
    pass


def rdkit_to_ase(rdG):
    block = Chem.MolToMolBlock(rdG)

    positions = np.empty((rdG.GetNumAtoms(), 3))
    symbols = []
    for i, atom in enumerate(block.split('\n')[4:rdG.GetNumAtoms() + 4]):
        data = atom.split()
        positions[i] = np.array(data[:3], dtype=float)
        symbols += [data[3]]

    return Atoms(symbols, positions)


def get_unsaturated_nodes(molecule, screen=None):

    unsaturated = []
    for node, data in molecule.nodes(data=True):
        radicals = data['radicals']

        if screen in data:
            continue

        if radicals > 0:
            unsaturated += [node]

    return np.array(unsaturated)


def isomorphic_molecules(molecule1, molecule2):
    em = iso.numerical_edge_match('bonds', 1)
    nm = iso.numerical_node_match('atomic_number', 1)
    isomorphic = nx.is_isomorphic(
        molecule1,
        molecule2,
        edge_match=em,
        node_match=nm)

    return isomorphic


def get_neighbor_symbols(G, u):
    """ Get chemical symbols of all neighboring atoms.
    """

    neighbor_symbols = []
    for v in G.neighbors(u):
        neighbor_symbols += [G.nodes(data=True)[v]['symbol']]
    return neighbor_symbols


def plot_molecule(molecule, file_name=None):

    G = get_rdkit_graph(molecule)
    G = Chem.RemoveHs(G)
    MolToFile(G, file_name, size=(200, 200))


def get_chemical_descriptor(graph, rank=1):

    symbols = nx.get_node_attributes(graph, 'symbol')
    chemical_tag = ''.join(np.sort(list(symbols.values())))

    if rank == 1:

        bonds = ''
        for u, v, data in graph.edges(data=True):
            bond = symbols[u] + symbols[v]
            bonds += bond * data.get('bonds')

        bonds_tag = ''.join(sorted(bonds))

        return chemical_tag, bonds_tag

    return chemical_tag


def get_rdkit_graph(molecule, sanitize=True):

    G = Chem.rdchem.EditableMol(Chem.rdchem.Mol())

    for j, data in molecule.nodes(data=True):
        rdAtom = Chem.rdchem.Atom(data['symbol'])
        rdAtom.SetNumRadicalElectrons(int(data['radicals']))
        G.AddAtom(rdAtom)

    rdBonds = Chem.rdchem.BondType
    orders = {
        '1': rdBonds.SINGLE,
        '2': rdBonds.DOUBLE,
        '3': rdBonds.TRIPLE}

    for u, v, data in molecule.edges(data=True):

        order = orders[str(data['bonds'])]
        G.AddBond(int(u), int(v), order)

    G = G.GetMol()

    if sanitize:
        Chem.SanitizeMol(G)

    return G


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

        Parameters:
            db_name (str): Name of the database file to access. Will
            connect to 'reaction-network.db' by default.
        """

        self.db_name = db_name
        self.molecules = {}
        self.base_free_radicals = {
            'H': 1,
            'C': 4,
            'O': 2,
            'N': 3,
            'S': 2}

        self.max_bond_limits = {
            'C': {
                'C': 3,
                'O': 2,
                'H': 1,
                'N': 3,
                'S': 2},
            'O': {
                'O': 2,
                'C': 2,
                'H': 1,
                'N': 2,
                'S': 2},
            'H': {
                'O': 1,
                'C': 1,
                'H': 1,
                'N': 1,
                'S': 1},
            'N': {
                'O': 2,
                'C': 2,
                'H': 1,
                'N': 2,
                'S': 2},
            'S': {
                'H': 1,
                'C': 2,
                'O': 2,
                'N': 2,
                'S': 2}}

        # Prevent formation of bond type YX
        # If X is bonded to species in list
        self.prevent_bond_chains = {}
        #     'OO': ['C', 'N', 'S'],
        #     'CO': ['O', 'N', 'S'],
        #     'SS': ['C', 'N', 'O'],
        #     'CS': ['S', 'N'],
        #     'NN': ['C', 'O', 'S'],
        #     'CN': ['N', 'S'],
        #     'SO': ['C', 'N', 'S', 'O'],
        #     'NO': ['C', 'N', 'S', 'O'],
        #     'OS': ['C', 'N', 'S', 'O'],
        #     'ON': ['C', 'N', 'S', 'O'],
        #     'SN': ['C', 'N', 'S', 'O'],
        #     'NS': ['C', 'N', 'S', 'O']
        # }

        # This will make terminations symmetric
        # tpbc = self.prevent_bond_chains.copy()
        # for YX, Z in tpbc.items():
        #     y, x = YX
        #     for z in Z:
        #         rchain = z + x
        #         if rchain not in self.prevent_bond_chains:
        #             self.prevent_bond_chains[rchain] = [y]
        #         elif y not in self.prevent_bond_chains[rchain]:
        #             self.prevent_bond_chains[rchain] += [y]

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
        This includes 4 tables: molecules, reactions, atoms, and bonds.
        """

        self.c.execute("""CREATE TABLE IF NOT EXISTS molecules(
        molecule_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        chemical_tag CHAR NOT NULL,
        bonds_tag CHAR NOT NULL
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS reactions(
        reaction_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        reactant1 INT NOT NULL,
        reactant2 INT NOT NULL,
        product1 INT NOT NULL,
        product2 INT NOT NULL,
        broken_bond CHAR,
        FOREIGN KEY(reactant1) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(reactant2) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(product1) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(product2) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(broken_bond) REFERENCES bonds(bond_pid),
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
        position_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INT NOT NULL,
        atom_id INT NOT NULL,
        x_coord REAL NOT NULL,
        y_coord REAL NOT NULL,
        z_coord REAL NOT NULL,
        symbol CHAR NOT NULL,
        FOREIGN KEY(molecule_id) REFERENCES molecules(molecule_pid),
        FOREIGN KEY(atom_id) REFERENCES atoms(atom_pid)
        )""")

        self.c.execute("""CREATE TABLE IF NOT EXISTS energies(
        energy_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        molecule_id INT NOT NULL,
        energy REAL NOT NULL,
        calculator CHAR,
        FOREIGN KEY(molecule_id) REFERENCES molecules(molecule_pid)
        )""")

    def molecule_search(
            self,
            element_pool={'C': 2, 'H': 6},
            load_molecules=True,
            multiple_bond_search=False):

        self.element_pool = element_pool
        self.multiple_bond_search = multiple_bond_search

        molecules = {}

        if load_molecules:
            self.molecules = self.load_molecules(binned=True)

        search_molecules = []
        for el in self.element_pool:

            molecule = nx.Graph()
            molecule.add_node(
                0,
                symbol=el,
                atomic_number=an[el],
                radicals=self.base_free_radicals[el])

            search_molecules += [molecule]

            if el not in self.molecules:
                molecules[el] = {'': [molecule]}

        for molecule in search_molecules:

            # Recusive use of molecules
            new_molecules, molecules = self.branch_molecule(
                molecule,
                molecules)

            search_molecules += new_molecules

        return molecules

    def branch_molecule(self, molecule, molecules):

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

                # Don't form new bonds if the users rules prohibit
                if self.max_bond_limits[base_el][el] == 0:
                    continue

                # Prevent user prohibited chains
                new_bond = ''.join(sorted(base_el + el))
                if new_bond in self.prevent_bond_chains:
                    neighbor_symbols = get_neighbor_symbols(
                        molecule,
                        base_node)

                    terminate = [
                        _ in neighbor_symbols for _ in
                        self.prevent_bond_chains[new_bond]]
                    if True in terminate:
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

                    if self.maximum_bond_limit(
                            molecule,
                            base_node,
                            existing_node):
                        continue

                    # Prevent user prohibited chains
                    new_bond = ''.join(sorted(base_el + el))
                    if new_bond in self.prevent_bond_chains:
                        neighbor_symbols = get_neighbor_symbols(
                            molecule,
                            existing_node)

                        term = [_ in neighbor_symbols for _ in
                                self.prevent_bond_chains[new_bond]]
                        if True in term:
                            continue

                    G = molecule.copy()

                    G.node[base_node]['radicals'] -= 1
                    G.node[existing_node]['radicals'] -= 1
                    if G.has_edge(base_node, existing_node):
                        G[base_node][existing_node]['bonds'] += 1
                    else:
                        G.add_edge(base_node, existing_node, bonds=1)

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
            reconfiguration=True,
            substitution=True):

        self.reconfiguration = reconfiguration

        molecules = self.load_molecules(binned=True)

        pathways, bbonds = [], []
        for data in molecules.values():
            for molecule_list in data.values():
                for molecule in molecule_list:

                    add_pathways, broken_bonds = self.get_addition_paths(
                        molecule,
                        molecules)

                    pathways += add_pathways
                    bbonds += broken_bonds

        self.save_pathways(pathways, bbonds)

        if reconfiguration:
            re_pathways = self.get_reconfiguration_paths(molecules, pathways)
            self.save_pathways(re_pathways)

        if substitution:
            sub_pathways = self.get_substitution_paths(molecules, pathways)
            self.save_pathways(sub_pathways)

        return pathways

    def get_addition_paths(
            self,
            molecule,
            molecules):

        disjoints, pathways, broken_bonds = [], [], []
        for u, v, data in molecule.edges(data=True):

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
                broken_bonds += ['{},{}'.format(u, v)]

                product_index = cut_molecule.graph['index']
                pieces = list(nx.connected_component_subgraphs(cut_molecule))

                addition_pathway = np.array([
                    [0, 0],
                    [product_index, 0]])
                for i, piece in enumerate(pieces):

                    chemical_tag, bonds_tag = get_chemical_descriptor(piece)

                    for reactant in molecules[chemical_tag][bonds_tag]:
                        if isomorphic_molecules(piece, reactant):
                            pindex = reactant.graph['index']
                            addition_pathway[0][i] = pindex
                            break

                [_.sort() for _ in addition_pathway]
                pathways += [addition_pathway]

        return pathways, broken_bonds

    def save_molecules(self, molecules):

        for chemical_tag, data in list(molecules.items()):
            for bonds_tag, molecule_list in list(data.items()):
                for molecule in molecule_list:

                    self.c.execute(
                        """INSERT INTO molecules
                        (chemical_tag, bonds_tag)
                        VALUES(?, ?)""", (chemical_tag, bonds_tag))

                    self.c.execute("""SELECT last_insert_rowid()""")
                    molecule_pid = self.c.fetchone()[0]

                    for u, v, data in molecule.edges(data=True):
                        bonds = data['bonds']
                        self.c.execute(
                            """INSERT INTO bonds
                            (molecule_id, node_id1, node_id2, nbonds)
                            VALUES(?, ?, ?, ?)""", (molecule_pid, u, v, bonds))

                    for node, data in molecule.nodes(data=True):
                        number = data['atomic_number']
                        radicals = data['radicals']
                        self.c.execute(
                            """INSERT INTO atoms
                            (molecule_id, node_id, atom_num, symbol, radicals)
                            VALUES(?, ?, ?, ?, ?)""",
                            (molecule_pid,
                             node,
                             number,
                             cs[number],
                             radicals))

    def save_pathways(self, pathways, broken_bonds=None):

        for i, path in enumerate(pathways):
            R, P = path
            P1, P2 = P
            R1, R2 = R

            if broken_bonds is not None:
                bbond = broken_bonds[i]
            else:
                bbond = None

            try:
                self.c.execute(
                    """INSERT INTO reactions
                    (product1, product2, reactant1, reactant2, broken_bond)
                    VALUES(?, ?, ?, ?, ?)""",
                    (int(P1), int(P2), int(R1), int(R2), bbond))
            except(sqlite3.IntegrityError):
                pass

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

    def load_pathways(self, broken_bonds=False):

        if broken_bonds:
            cmd = """SELECT
            reactant1,
            reactant2,
            product1,
            product2,
            broken_bond
            FROM reactions"""
        else:
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
            _add = []
            for _ in path:
                try:
                    _add += [int(_)]
                except(ValueError, TypeError):
                    _add += [_]

            pathways += [_add]

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
                    if R and P:
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

    def maximum_bond_limit(self, molecule, n1, n2):

        el1 = molecule.node[n1]['symbol']
        el2 = molecule.node[n2]['symbol']

        bonds = 0
        if molecule.has_edge(n1, n2):
            bonds = molecule[n1][n2]['bonds']

        if self.max_bond_limits[el1][el2] == bonds:
            return True

        return False

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
            if cycles > 0:
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

        for u, v, data in molecule.edges(data=True):
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

        new_pathways = []
        reconfig = set()
        for R, P in pathways:
            R1, R2 = R
            P1, P2 = P

            if not R1 or P1:
                continue

            mol_R1 = ind_mol[R1]
            mol_R2 = ind_mol[R2]

            reconfigurations = [ind_mol[P2]]

            # Find potential bonding sites
            p1_bonds = get_unsaturated_nodes(
                ind_mol[R1])

            p2_bonds = get_unsaturated_nodes(
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
            reconfig_pathway = np.array([
                [P2, 0],
                [0, 0]])
            for Pt in reconfigurations:
                chemical_tag, bonds_tag = get_chemical_descriptor(Pt)

                # Bond structure has not been enumerated
                if bonds_tag not in molecules[chemical_tag]:
                    continue

                for P in molecules[chemical_tag][bonds_tag]:
                    if isomorphic_molecules(P, Pt):
                        reconfig_pathway[1][0] = P.graph['index']
                        pc = tuple(sorted([P.graph['index'], P2]))
                        break

                if pc not in reconfig:
                    reconfig.add(pc)
                    new_pathways += [reconfig_pathway]

        for path in new_pathways:
            [path.sort() for path in new_pathways]
        return new_pathways

    def save_3d_structure(self, molecule, uff=0):
        G = get_rdkit_graph(molecule)
        Chem.EmbedMolecule(G, Chem.ETKDG())

        lec = 0
        if uff:
            cids = Chem.EmbedMultipleConfs(G, numConfs=uff)
            Chem.UFFOptimizeMoleculeConfs(G)

            energies = []
            for cid in cids:
                ffe = Chem.UFFGetMoleculeForceField(G, confId=cid).CalcEnergy()
                energies += [ffe]
            energies = np.array(energies)

            lec = np.argmin(energies)

        block = Chem.MolToMolBlock(G, confId=int(lec))
        for j, atom in enumerate(block.split('\n')[4:G.GetNumAtoms()+4]):
            data = atom.split()
            x, y, z = np.array(data[:3], dtype=float)

            self.c.execute(
                """INSERT INTO positions
                (molecule_id, atom_id, x_coord, y_coord, z_coord, symbol)
                VALUES(?, ?, ?, ?, ?, ?)""",
                (molecule.graph['index'],
                 j, x, y, z, data[3]))

    def load_3d_structures(self, molecule_id=None):
        images = []
        if isinstance(molecule_id, list):
            molecule_id = ','.join([str(_) for _ in molecule_id])

        if molecule_id is None:
            cmd = """SELECT
             GROUP_CONCAT(x_coord || ',' || y_coord || ',' || z_coord, ';'),
             GROUP_CONCAT(symbol, ';')
             FROM positions
             GROUP BY molecule_id
            """
        else:
            cmd = """SELECT
             GROUP_CONCAT(x_coord || ',' || y_coord || ',' || z_coord, ';'),
             GROUP_CONCAT(symbol, ';')
             FROM positions
             WHERE molecule_id IN ({})
             GROUP BY molecule_id
            """.format(molecule_id)

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        for out in fetch:

            symbols = out[1].split(';')
            positions = np.array(
                [_.split(',')
                 for _ in out[0].split(';')], dtype=float)

            images += [Atoms(symbols, positions)]

        return images

    def get_substitution_paths(
            self,
            molecules,
            pathways):
        # Follows the form:
        # R1(-P1) + R2(-P1) --> P1 + P2

        substitution_set = set()
        new_pathways = []

        ind_mol = self.load_molecules()

        # Get the maximum number of nodes
        maxn = 0
        for i in ind_mol:
            G = ind_mol[i]
            n = G.number_of_nodes()
            if n > maxn:
                maxn = n

        # Need to eliminate isomorphs to
        # not count double bonds
        for iP1 in ind_mol:

            P1 = ind_mol[iP1]
            P1_bonds = get_unsaturated_nodes(P1)
            nP1 = P1.number_of_nodes()

            if len(P1_bonds) == 0:
                continue

            for iRa, iPa in pathways:
                iR1, iR2 = iRa
                _, iP2 = iPa

                # Ignore bonding pathways
                if not iR1:
                    continue

                for i, iR in enumerate(iRa):
                    R = ind_mol[iR]
                    nR = R.number_of_nodes()

                    # Don't create larger molecules
                    if nR + nP1 > maxn:
                        continue

                    R_bonds = get_unsaturated_nodes(R) + nP1

                    for b1 in P1_bonds:
                        for b2 in R_bonds:

                            R_P1 = nx.disjoint_union(P1, R)
                            R_P1.add_edge(b1, b2, bonds=1)

                            (chemical_tag,
                             bonds_tag) = get_chemical_descriptor(R_P1)

                            # Bond structure has not been enumerated
                            if chemical_tag not in molecules:
                                continue

                            # Bond structure has not been enumerated
                            if bonds_tag not in molecules[chemical_tag]:
                                continue

                            subst_pathway = np.array([
                                sorted([0, iRa[::-1][i]]),
                                sorted([iP1, iP2])])

                            for G in molecules[chemical_tag][bonds_tag]:
                                if isomorphic_molecules(R_P1, G):
                                    iR_P1 = G.graph['index']

                                    subst_pathway = np.array([
                                        sorted([iR_P1, iRa[::-1][i]]),
                                        sorted([iP1, iP2])])

                                    sum_index = (
                                        ','.join(subst_pathway[0].astype(str)),
                                        ','.join(subst_pathway[1].astype(str)))

                                    if sum_index not in substitution_set:
                                        # These are identical
                                        substitution_set.add(sum_index)
                                        substitution_set.add(sum_index[::-1])

                                        new_pathways += [subst_pathway]
                                    break

        return new_pathways
