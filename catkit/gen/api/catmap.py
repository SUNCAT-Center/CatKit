from ..pathways import ReactionNetwork
from .rd_kit import get_smiles
import numpy as np
import networkx as nx
from ase.utils import formula_hill
from collections import Counter


class catmap_api():
    """This class contains information and functions for exporting
    input to catmap."""

    def __init__(self,
                 pathways=None,
                 molecules=None,
                 db_name=None,
                 formula='hill'):
        if pathways is not None and molecules is not None:
            self.pathways = pathways
            self.molecules = molecules
        elif db_name is not None:
            with ReactionNetwork(db_name=db_name) as rn:
                self.pathways = rn.path_search(
                    reconfiguration=True, substitution=True)
                self.molecules = rn.load_molecules()
        else:
            raise ValueError("Pass a list of reaction edges in 'pathways'" +
                             " and a list of molecule graphs in 'molecules'," +
                             " or the path and filename to a sqlite3 file in" +
                             " in 'db_name'.")
        self.formula = formula

    def rxn_expressions(self):
        """Returns a list of reaction expressions in catmap format."""
        rxn_expressions = []
        # Iterate over reactions and remove placeholders.
        for r in self.pathways:
            R = r[0]
            R = R[np.nonzero(R)]
            P = r[1]
            P = P[np.nonzero(P)]
            # Append individual reaction strings to a list.
            rxn_expressions.append(self._rxn_equation(R, P))
        return rxn_expressions

    def species_definitions(self):
        species_definitions = {}
        for index in self.molecules:
            species_definitions.update(self._species_definition(index))
        return species_definitions

    def _rxn_equation(self, R, P, ts=False, formula='hill'):
        """Return a reaction equation in catmap format.
        Parameters
        ----------
        R : 2 x 2 array
            2 edges identifying the initial state
        P : 2 x 2 array
            2 edges identifying the final state
        """
        # Ensure the reaction is balanced.
        R, P = self._balance_equation(R, P)
        reactants = self._get_state(R)
        products = self._get_state(P)
        if ts:
            raise NotImplementedError("Transition states.")
            # return (' + '.join(reactants) +
            #         ' <-> ' + ' + '.join(transition_state) +
            #         ' -> ' + '.join(products))
        else:
            return ' + '.join(reactants) + ' <-> ' + ' + '.join(products)

    def _species_definition(self, index):
        """Return a dictionary with species definitions for catmap.
        Catmap needs information about the composition of species,
        if they are not named in either hill notation or
        condensed structural notation."""
        symbols = nx.get_node_attributes(self.molecules[index], 'symbol')
        c = Counter({})
        for s in symbols.values():
            c += Counter({s: 1})
        return {self._get_species_name(index): {'composition': dict(c)}}

    def _get_state(self, index):
        """Return a list of species names from a list of molecular indices.
        Parameters
        ----------
        index : list
            list of integers referring to elements in self.molecules.
        formula : str
            'hill' or 'smiles'
        """
        state = []
        for mol in index:
            species = self._get_species_name(mol)
            state.append(species)
        return state

    def _get_species_name(self, index):
        site = self._get_site_index(index)
        if self.formula == 'hill':
            n = nx.get_node_attributes(self.molecules[index], 'atomic_number')
            name = formula_hill(list(n.values()))
        elif self.formula == 'smiles':
            name = get_smiles(self.molecules[index])
        else:
            raise NotImplementedError(str(self.formula))
        return '_'.join([name, site])

    def _get_site_index(self, index):
        """Return the site of a molecule.
        Parameters
        ----------
        molecule : networkx graph object
            CatKit molecule.
        Todo:
            implement surface sites in CatKit molecules.
        """
        return 'g'

    def _balance_equation(self, R, P):
        """Returns a balanced reaction equation."""
        c_R = self._count_atoms(R)
        c_P = self._count_atoms(P)
        if c_R == c_P:
            return R, P
        else:
            print(c_R, c_P)
            raise NotImplementedError("Return balanced equation.")

    def _count_atoms(self, index):
        """Count the atoms in a given state."""
        c = Counter({})
        for mol in index:
            symbols = nx.get_node_attributes(self.molecules[mol], 'symbol')
            for s in symbols.values():
                c += Counter({s: 1})
        return c
