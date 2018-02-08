import sqlite3
import numpy as np
import networkx as nx
from ase import Atom
import networkx.algorithms.isomorphism as iso
from ase.data import chemical_symbols as cs
import matplotlib.pyplot as plt
from catkit import Gratoms
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

    return Gratoms(symbols, positions)


def plot_molecule(molecule, file_name=None):

    G = get_rdkit_graph(molecule)
    G = Chem.RemoveHs(G)
    MolToFile(G, file_name, size=(200, 200))


def get_rdkit_graph(molecule, sanitize=True):

    G = Chem.rdchem.EditableMol(Chem.rdchem.Mol())

    for j, data in molecule.nodes(data=True):
        rdAtom = Chem.rdchem.Atom(cs[data['number']])
        rdAtom.SetNumRadicalElectrons(int(data['valence']))
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
    """ """
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
            db_name='reaction-network.db',
            base_valence=None,
            nbond_limits=None,
    ):
        """ The __init__ function is automatically called when the
        class is referenced.

        Parameters:
            db_name (str): Name of the database file to access. Will
            connect to 'reaction-network.db' by default.
        """

        self.db_name = db_name
        self.molecules = {}

        if base_valence is None:
            base_valence = np.zeros(17)
            base_valence[[1, 6, 7, 8, 16]] = [1, 4, 3, 2, 2]
            self.base_valence = base_valence

        if nbond_limits is None:
            nbond_limits = np.ones((17, 17))
            nbond_limits[1:, 1:] = 2
            nbond_limits[6, [6, 7]] = 3
            self.nbond_limits = nbond_limits

    def __enter__(self):
        """ This function is automatically called whenever the class
        is used together with a 'with' statement.
        """
        self.con = sqlite3.connect(self.db_name)
        self.c = self.con.cursor()
        self.create_table()

        return self

    def __exit__(self, type, value, tb):
        """ Upon exiting the 'with' statement, __exit__ is called."""
        self.con.commit()
        self.con.close()

    def create_table(self):
        """ Creates the database table framework used in SQLite.
        This includes 4 tables: molecules, reactions, atoms, and bonds.
        """
        self.c.execute("""CREATE TABLE IF NOT EXISTS molecules(
        molecule_pid INTEGER PRIMARY KEY AUTOINCREMENT,
        comp_tag CHAR NOT NULL,
        bond_tag CHAR NOT NULL
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
        valence INT NOT NULL,
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
            multiple_bond_search=False
    ):

        numbers = np.zeros(17)
        for k, v in element_pool.items():
            numbers[cs.index(k)] = v

        self.element_pool = numbers
        self.multiple_bond_search = multiple_bond_search

        molecules = {}
        if load_molecules:
            self.molecules = self.load_molecules(binned=True)

        search_molecules = []
        for el, cnt in enumerate(self.element_pool):

            if cnt == 0:
                continue

            molecule = Gratoms(symbols=[el])
            molecule.nodes[0]['valence'] = self.base_valence[el]

            search_molecules += [molecule]

            comp_tag, bond_tag = molecule.get_chemical_tags()
            if comp_tag not in self.molecules:
                molecules[comp_tag] = {bond_tag: [molecule]}

        for molecule in search_molecules:

            # Recusive use of molecules
            new_molecules, molecules = self.branch_molecule(
                molecule,
                molecules
            )

            search_molecules += new_molecules

        return molecules

    def branch_molecule(self, molecule, molecules):

        new_molecules = []
        nodes = molecule.nodes
        counter = np.bincount(
            molecule.get_atomic_numbers(),
            minlength=17
        )

        # Creating new nodes
        for el, cnt in enumerate(self.element_pool):

            # Skip elements exceed specified limit
            if counter[el] >= cnt:
                continue

            for base_node in nodes:

                # Check if an additional bond can be formed
                valence = nodes[base_node]['valence']
                base_el = molecule.nodes[base_node]['number']

                if valence <= 0:
                    continue

                # Don't form new bonds if the users rules prohibit
                if self.nbond_limits[base_el, el] == 0:
                    continue

                G = molecule.copy()
                node_index = len(G)

                G.nodes[base_node]['valence'] -= 1
                G += Atom(el)
                G.nodes[node_index]['valence'] = self.base_valence[el] - 1
                G.graph.add_edge(base_node, node_index, bonds=1)

                comp_tag, bond_tag = G.get_chemical_tags()

                isomorph_found = False
                if comp_tag not in molecules:
                    molecules[comp_tag] = {}

                if bond_tag not in molecules[comp_tag]:
                    molecules[comp_tag][bond_tag] = []

                else:
                    for G1 in molecules[comp_tag][bond_tag]:
                        if G.is_isomorph(G1):
                            isomorph_found = True
                            break

                if not isomorph_found:
                    new_molecules += [G]
                    molecules[comp_tag][bond_tag] += [G]

            # Bonds to existing nodes
            if self.multiple_bond_search:
                for existing_node in nodes:

                    if valence <= 0:
                        continue

                    # Atoms don't bond to themselves
                    if existing_node == base_node:
                        continue

                    valence_new = nodes[existing_node]['valence']
                    el = molecule.nodes[existing_node]['number']

                    if valence_new <= 0:
                        continue

                    if self.maximum_bond_limit(
                            molecule,
                            base_node,
                            existing_node):
                        continue

                    G = molecule.copy()

                    G.nodes[base_node]['valence'] -= 1
                    G.nodes[existing_node]['valence'] -= 1
                    if G.graph.has_edge(base_node, existing_node):
                        G.graph[base_node][existing_node]['bonds'] += 1
                    else:
                        G.graph.add_edge(base_node, existing_node, bonds=1)

                    comp_tag, bond_tag = G.get_chemical_tags()

                    isomorph_found = False
                    if comp_tag not in molecules:
                        molecules[comp_tag] = {}

                    if bond_tag not in molecules[comp_tag]:
                        molecules[comp_tag][bond_tag] = []

                    else:
                        for G1 in molecules[comp_tag][bond_tag]:
                            if G.is_isomorph(G1):
                                isomorph_found = True
                                break

                    if not isomorph_found:
                        new_molecules += [G]
                        molecules[comp_tag][bond_tag] += [G]

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

            cut_molecule.graph[u][v]['bonds'] -= 1
            if not cut_molecule.graph[u][v]['bonds']:
                cut_molecule.graph.remove_edge(u, v)

            isomorph_found = False
            for disjoint in disjoints:
                if cut_molecule.is_isomorph(disjoint):
                    isomorph_found = True
                    break

            if not isomorph_found:
                disjoints += [cut_molecule]
                broken_bonds += ['{},{}'.format(u, v)]

                product_index = cut_molecule.graph.name
                pieces = list(nx.connected_components(cut_molecule.graph))

                addition_pathway = np.array([
                    [0, 0],
                    [product_index, 0]])
                for i, piece in enumerate(pieces):

                    if len(pieces) == 2:
                        G = molecule.copy()
                        del G[list(piece)]
                    else:
                        G = cut_molecule.copy()

                    comp_tag, bond_tag = G.get_chemical_tags()

                    for reactant in molecules[comp_tag][bond_tag]:
                        if G.is_isomorph(reactant):
                            pindex = reactant.graph.name
                            addition_pathway[0][i] = pindex
                            break

                [_.sort() for _ in addition_pathway]
                pathways += [addition_pathway]

        return pathways, broken_bonds

    def save_molecules(self, molecules):

        for comp_tag, data in list(molecules.items()):
            for bond_tag, molecule_list in list(data.items()):
                for molecule in molecule_list:

                    self.c.execute(
                        """INSERT INTO molecules
                        (comp_tag, bond_tag)
                        VALUES(?, ?)""", (comp_tag, bond_tag))

                    self.c.execute("""SELECT last_insert_rowid()""")
                    molecule_pid = self.c.fetchone()[0]

                    for u, v, data in molecule.edges(data=True):
                        bonds = data['bonds']
                        self.c.execute(
                            """INSERT INTO bonds
                            (molecule_id, node_id1, node_id2, nbonds)
                            VALUES(?, ?, ?, ?)""", (molecule_pid, u, v, bonds))

                    for node, data in molecule.nodes(data=True):
                        number = data['number']
                        valence = data['valence']
                        self.c.execute(
                            """INSERT INTO atoms
                            (molecule_id, node_id, atom_num, valence)
                            VALUES(?, ?, ?, ?)""",
                            (molecule_pid,
                             node,
                             int(number),
                             valence))

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
                        m.comp_tag,
                        m.bond_tag,
                        nodes,
                        bonds
        FROM molecules m LEFT JOIN
        (
         SELECT molecule_id,
         GROUP_CONCAT(node_id || ',' || atom_num || ',' || valence, ';')
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

        for index, comp_tag, bond_tag, node_data, edge_data in fetch:

            # Unpacks node, number, and valence
            node_data = np.array(
                [_.split(',') for _ in node_data.split(';')],
                dtype=int
            )

            data, symbols = {}, []
            for node, n, valence in node_data:
                data.update({node: valence})
                symbols += [n]

            molecule = Gratoms(symbols)
            molecule.graph.name = index
            nx.set_node_attributes(molecule.graph, data, 'valence')

            if edge_data:
                edges = np.array(
                    [_.split(',') for _ in edge_data.split(';')],
                    dtype=int
                )
                molecule.graph.add_weighted_edges_from(edges, weight='bonds')

            if binned:
                if comp_tag not in molecules:
                    molecules[comp_tag] = {}

                if bond_tag not in molecules[comp_tag]:
                    molecules[comp_tag][bond_tag] = []

                molecules[comp_tag][bond_tag] += [molecule]
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
            p1_bonds = ind_mol[R1].get_unsaturated_nodes()

            p2_bonds = ind_mol[R2].get_unsaturated_nodes() + \
                       len(ind_mol[R1])

            for b1 in p1_bonds:
                for b2 in p2_bonds:

                    Pt = mol_R1 + mol_R2
                    Pt.graph.add_edge(b1, b2, bonds=1)

                    isomorph_found = False
                    for P in reconfigurations:
                        if Pt.is_isomorph(P):
                            isomorph_found = True
                            break

                    if not isomorph_found:
                        reconfigurations += [Pt]

            del reconfigurations[0]
            reconfig_pathway = np.array([
                [P2, 0],
                [0, 0]])
            for Pt in reconfigurations:
                comp_tag, bond_tag = Pt.get_chemical_tags()

                # Bond structure has not been enumerated
                if bond_tag not in molecules[comp_tag]:
                    continue

                for P in molecules[comp_tag][bond_tag]:
                    if P.is_isomorph(Pt):
                        reconfig_pathway[1][0] = P.graph.name
                        pc = tuple(sorted([P.graph.name, P2]))
                        break

                if pc not in reconfig:
                    reconfig.add(pc)
                    new_pathways += [reconfig_pathway]

        for path in new_pathways:
            [path.sort() for path in new_pathways]
        return new_pathways

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
            n = len(G)
            if n > maxn:
                maxn = n

        # Need to eliminate isomorphs to
        # not count double bonds
        for iP1 in ind_mol:

            P1 = ind_mol[iP1]
            P1_bonds = P1.get_unsaturated_nodes()
            nP1 = len(P1)

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
                    nR = len(R)

                    # Don't create larger molecules
                    if nR + nP1 > maxn:
                        continue

                    R_bonds = R.get_unsaturated_nodes() + nP1

                    for b1 in P1_bonds:
                        for b2 in R_bonds:

                            R_P1 = P1 + R
                            R_P1.graph.add_edge(b1, b2, bonds=1)

                            comp_tag, bond_tag = R_P1.get_chemical_tags()

                            # Bond structure has not been enumerated
                            if comp_tag not in molecules or \
                               bond_tag not in molecules[comp_tag]:
                                continue

                            subst_pathway = np.array([
                                sorted([0, iRa[::-1][i]]),
                                sorted([iP1, iP2])])

                            for G in molecules[comp_tag][bond_tag]:
                                if R_P1.is_isomorph(G):
                                    iR_P1 = G.graph.name

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
                (molecule.graph.name,
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

            images += [Gratoms(symbols, positions)]

        return images

    def plot_reaction_network(self, file_name=None):

        pathways = []
        load_paths = self.load_pathways()
        for path in load_paths:
            for R in path[:2]:
                for P in path[2:]:
                    if R and P:
                        pathways += [[R, P]]

        network = nx.Graph(pathways)

        plt.figure(figsize=(6, 6))
        nx.draw_networkx(network)
        plt.draw()
        plt.axis('off')

        if file_name:
            plt.savefig(file_name)

        plt.close()

    def maximum_bond_limit(self, molecule, n1, n2):

        el1 = molecule.nodes[n1]['number']
        el2 = molecule.nodes[n2]['number']

        bonds = 0
        if molecule.graph.has_edge(n1, n2):
            bonds = molecule.graph[n1][n2]['bonds']

        if self.nbond_limits[el1][el2] == bonds:
            return True

        return False
