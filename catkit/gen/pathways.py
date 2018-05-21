from catkit import Gratoms
import sqlite3
import numpy as np
import networkx as nx
from ase import Atom
from ase.data import chemical_symbols
import matplotlib.pyplot as plt


class ReactionNetwork():
    """ A class for accessing a temporary SQLite database. This
    function works as a context manager and should be used as follows:

    with ReactionNetwork() as rn:
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
        """The __init__ function is automatically called when the
        class is referenced.

        Parameters
        ----------
        db_name : str
            Name of the database file to access. Will connect to
            'reaction-network.db' by default.
        base_valence : array (n,)
            The maximum number of bonds each chemical species may form. Each
            entry of the array is the bond limit of the corresponding atomic
            number. The 0th entry of the array is ignored.
        nbond_limits : array (n, n):
            The maximum number of bonds that can be formed when species of
            index 0 is bonding to a species of index 1. Each entry of the array
            is the bond limit of the corresponding atomic number. The 0th entry
            of the array is ignored.
        """

        self.db_name = db_name
        self.molecules = {}

        if base_valence is None:
            base_valence = np.zeros(17)
            # Define valence for H, C, N, O, and S
            base_valence[[1, 6, 7, 8, 16]] = [1, 4, 3, 2, 2]
            self.base_valence = base_valence

        if nbond_limits is None:
            nbond_limits = np.ones((17, 17))
            nbond_limits[2:, 2:] = 2
            # # Specific bonding limits for Carbon
            nbond_limits[6, [6, 7]] = 3
            self.nbond_limits = nbond_limits

    def __enter__(self):
        """Initialize the database connection when entering
        an indented loop.
        """
        self.con = sqlite3.connect(self.db_name)
        self.c = self.con.cursor()
        self.create_table()

        return self

    def __exit__(self, type, value, tb):
        """Commit and close the database upon exiting indentation."""
        self.con.commit()
        self.con.close()

    def create_table(self):
        """Create the SQLite database table framework."""

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

    def molecule_search(self,
                        element_pool={'C': 2, 'H': 6},
                        load_molecules=True,
                        multiple_bond_search=False):
        """Return the enumeration of molecules which can be produced from
        the specified atoms.

        Parameters
        ----------
        element_pool : dict
            Atomic symbols keys paired with the maximum number of that atom.
        load_molecules : bool
            Load any existing molecules from the database.
        multiple_bond_search : bool
            Allow atoms to form bonds with other atoms in the molecule.
        """
        numbers = np.zeros(len(self.base_valence))
        for k, v in element_pool.items():
            numbers[chemical_symbols.index(k)] = v

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
            new_molecules, molecules = self._branch_molecule(
                molecule, molecules)

            search_molecules += new_molecules

        self.save_molecules(molecules)

    def _branch_molecule(self, molecule, molecules):
        """Construct all the molecules one atom larger from a given
        molecule and add the unique ones to a dictionary.

        Parameters:
        -----------
        molecule : Gratoms object
            Molecule to create branching structures from.
        molecules : dict
            Molecule structures to check for unique matches against.

        Returns:
        --------
        new_molecules : list
            All unique molecules discovered while branching.
        molecules : dict
            Molecule structures to check for unique matches updated with
            new_molecules.
        """
        new_molecules = []
        nodes = molecule.nodes
        counter = np.bincount(
            molecule.get_atomic_numbers(), minlength=len(self.base_valence))

        for base_node in nodes:

            # Check if an additional bond can be formed
            valence = nodes[base_node]['valence']
            base_el = molecule.nodes[base_node]['number']

            # Creating new nodes
            for el, cnt in enumerate(self.element_pool):

                # Skip elements exceeding specified limit
                if counter[el] >= cnt:
                    continue

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

                    if self._maximum_bond_limit(molecule, base_node,
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

    def path_search(self, reconfiguration=True, substitution=True):
        """Search for reaction mechanisms from a pre-populated database
        of molecules. By default, only single bond addition pathways are
        enumerated (Also called elementary steps).

        Parameters
        ----------
        reconfiguration : bool
            Search for reconfiguration paths. Reconfiguration paths are
            all those where only the bond order is changed. R1 --> P1.
        substitution : bool
            Search for substitution paths. Substitution paths are all
            those where one bond is broken and one bond is formed
            simultaneously. R1 + R2 --> P1 + P2.
        """
        molecules = self.load_molecules(binned=True)

        pathways, bbonds = [], []
        for data in molecules.values():
            for molecule_list in data.values():
                for molecule in molecule_list:

                    add_pathways, broken_bonds = self._get_addition_paths(
                        molecule, molecules)

                    pathways += add_pathways
                    bbonds += broken_bonds

        self.save_pathways(pathways, bbonds)

        if reconfiguration:
            re_pathways = self._get_reconfiguration_paths(molecules, pathways)
            self.save_pathways(re_pathways)

        if substitution:
            sub_pathways = self._get_substitution_paths(molecules, pathways)
            self.save_pathways(sub_pathways)

    def _get_addition_paths(self, molecule, molecules):
        """TODO: I can haz documentation?"""
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

                addition_pathway = np.array([[0, 0], [product_index, 0]])
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

    def _get_reconfiguration_paths(self, molecules, pathways):
        """TODO: I can haz documentation?"""
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
            reconfig_pathway = np.array([[P2, 0], [0, 0]])
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

    def _get_substitution_paths(self, molecules, pathways):
        """TODO: I can haz documentation?

        Follows the form:
        R1(-P1) + R2(-P1) --> P1 + P2
        """
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
                                sorted([iP1, iP2])
                            ])

                            for G in molecules[comp_tag][bond_tag]:
                                if R_P1.is_isomorph(G):
                                    iR_P1 = G.graph.name

                                    subst_pathway = np.array([
                                        sorted([iR_P1, iRa[::-1][i]]),
                                        sorted([iP1, iP2])
                                    ])

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

    def save_molecules(self, molecules):
        """Save enumerated molecules to the ReactionNetwork database.

        Parameters
        ----------
        molecules : dict
            Molecules to be saved to the database.
        """
        for comp_tag, data in molecules.items():
            for bond_tag, molecule_list in data.items():
                for molecule in molecule_list:

                    self.c.execute("""INSERT INTO molecules
                        (comp_tag, bond_tag)
                        VALUES(?, ?)""", (comp_tag, bond_tag))

                    self.c.execute("""SELECT last_insert_rowid()""")
                    molecule_pid = self.c.fetchone()[0]

                    for u, v, data in molecule.edges(data=True):
                        bonds = data['bonds']
                        self.c.execute(
                            """INSERT INTO bonds
                            (molecule_id, node_id1, node_id2, nbonds)
                            VALUES(?, ?, ?, ?)""",
                            (molecule_pid, int(u), int(v), int(bonds)))

                    for node, data in molecule.nodes(data=True):
                        number = int(data.get('number'))
                        valence = data.get('valence')

                        if valence is None:
                            valence = -1

                        self.c.execute(
                            """INSERT INTO atoms
                            (molecule_id, node_id, atom_num, valence)
                            VALUES(?, ?, ?, ?)""",
                            (molecule_pid, int(node), number, int(valence)))

    def save_pathways(self, pathways, broken_bonds=None):
        """Save enumerated pathways the ReactionNetwork database.
        More than two reactants or two products is not supported.

        Parameters
        ----------
        pathways : list
            Sorted pathways in the form [R1, R2, P1, P2].
        broken_bonds : list
            Comma separated strings of index associated with the
            two atoms whos bond is broken. List order must match
            pathways.
        """
        for i, path in enumerate(pathways):
            R, P = path
            P1, P2 = P
            R1, R2 = R

            if broken_bonds is not None:
                bbond = broken_bonds[i]
            else:
                bbond = None

            try:
                self.c.execute("""INSERT INTO reactions
                    (product1, product2, reactant1, reactant2, broken_bond)
                    VALUES(?, ?, ?, ?, ?)""",
                               (int(P1), int(P2), int(R1), int(R2), bbond))
            except (sqlite3.IntegrityError):
                pass

    def load_molecules(self, ids=None, binned=False):
        """Load 2D molecule graphs from the database.

        Parameters
        ----------
        binned : bool
            Return the molecules in sub-dictionaries of their
            corresponding composition and bonding tags.

        Returns
        -------
        molecules : dict
            All molecules present in the database.
        """
        if isinstance(ids, list):
            ids = ','.join([str(_) for _ in ids])

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
        if ids:
            cmd += """WHERE m.molecule_pid IN ({})""".format(ids)

        self.c.execute(cmd)
        fetch = self.c.fetchall()

        molecules = {}

        for index, comp_tag, bond_tag, node_data, edge_data in fetch:

            # Unpacks node, number, and valence
            node_data = np.array(
                [_.split(',') for _ in node_data.split(';')], dtype=int)

            data, symbols = {}, []
            for node, n, valence in node_data:
                data.update({node: valence})
                symbols += [n]

            molecule = Gratoms(symbols)
            molecule.graph.name = index
            nx.set_node_attributes(molecule.graph, data, 'valence')

            if edge_data:
                edges = np.array(
                    [_.split(',') for _ in edge_data.split(';')], dtype=int)
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
        """Save enumerated pathways the ReactionNetwork database.
        More than two reactants or two products is not supported.

        Parameters
        ----------
        broken_bonds : bool
            Return the index information of which bond was broken.
            Only supported for elementary steps.

        Returns
        -------
        pathways : list
            All pathways present in the database.
        """
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
                except (ValueError, TypeError):
                    _add += [_]

            pathways += [_add]

        pathways = np.array(pathways)

        return pathways

    def save_3d_structure(self, gratoms, overwrite=False):
        """Save Cartesian coordinates into the ReactionNetwork database.

        Parameters
        ----------
        gratoms : Gratoms or Atoms object
            Object with Cartesian coordinates to be saved.
        overwrite : bool
            Allow the database to overwrite a matching index.
        """
        name = gratoms.graph.name

        if overwrite:
            cmd = """SELECT GROUP_CONCAT(position_pid)
            FROM positions
            WHERE molecule_id = ({})
            GROUP BY molecule_id
            """.format(name)

            self.c.execute(cmd)
            match = self.c.fetchone()

            if match:
                match = match[0].split(',')

                assert (len(match) == len(gratoms))

                for i, atom in enumerate(gratoms):
                    x, y, z = atom.position

                    self.c.execute("""UPDATE positions
                        SET x_coord = ?,
                            y_coord = ?,
                            z_coord = ?
                        WHERE position_pid = ?
                        """, (x, y, z, match[i]))

                return

        for i, atom in enumerate(gratoms):
            x, y, z = atom.position
            symbol = atom.get('symbol')

            self.c.execute("""INSERT INTO positions
                (molecule_id, atom_id, x_coord, y_coord, z_coord, symbol)
                VALUES(?, ?, ?, ?, ?, ?)""", (name, i, x, y, z, symbol))

    def load_3d_structures(self, ids=None):
        """Return Gratoms objects from the ReactionNetwork database.

        Parameters
        ----------
        ids : int or list of int
            Identifier of the molecule in the database. If None, return all
            structure.

        Returns
        -------
        images : list
            All Gratoms objects in the database.
        """
        if isinstance(ids, list):
            ids = ','.join([str(_) for _ in ids])

        if ids is None:
            cmd = """SELECT
             GROUP_CONCAT(x_coord || ',' || y_coord || ',' || z_coord, ';'),
             GROUP_CONCAT(symbol, ';')
             FROM positions
             GROUP BY molecule_id
            """

            self.c.execute(cmd)
            fetch = self.c.fetchall()

            images = []
            for i, out in enumerate(fetch):

                symbols = out[1].split(';')
                positions = np.array(
                    [_.split(',') for _ in out[0].split(';')], dtype=float)

                gratoms = Gratoms(symbols, positions)
                gratoms.graph.name = i

                images += [gratoms]

            return images

        else:
            cmd = """SELECT
             GROUP_CONCAT(x_coord || ',' || y_coord || ',' || z_coord, ';'),
             GROUP_CONCAT(symbol, ';')
             FROM positions
             WHERE molecule_id IN ({})
             GROUP BY molecule_id
            """.format(ids)

            self.c.execute(cmd)
            out = self.c.fetchone()

            if out is None:
                raise ValueError('No matching index found')

            symbols = out[1].split(';')
            positions = np.array(
                [_.split(',') for _ in out[0].split(';')], dtype=float)

            gratoms = Gratoms(symbols, positions)
            gratoms.graph.name = int(ids)

            return gratoms

    def plot_reaction_network(self, file_name=None):
        """Plot the reaction network present in the database."""

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

    def _maximum_bond_limit(self, molecule, n1, n2):
        """Return whether a maximum bonding limit has been met."""
        el1 = molecule.nodes[n1]['number']
        el2 = molecule.nodes[n2]['number']

        bonds = 0
        if molecule.graph.has_edge(n1, n2):
            bonds = molecule.graph[n1][n2]['bonds']

        if self.nbond_limits[el1][el2] == bonds:
            return True

        return False
