from .. import Gratoms
from ase.data import chemical_symbols
import numpy as np
import networkx as nx
try:
    from rdkit.Chem.Draw import MolToFile
    from rdkit.Chem import AllChem as Chem
except(ImportError):
    pass


def get_graph(molecule, sanitize=True):

    rdkG = Chem.rdchem.EditableMol(Chem.rdchem.Mol())

    for j, data in molecule.nodes(data=True):
        rdAtom = Chem.rdchem.Atom(chemical_symbols[data['number']])
        rdAtom.SetNumRadicalElectrons(int(data['valence']))
        rdkG.AddAtom(rdAtom)

    rdBonds = Chem.rdchem.BondType
    orders = {
        '1': rdBonds.SINGLE,
        '2': rdBonds.DOUBLE,
        '3': rdBonds.TRIPLE
    }

    for u, v, data in molecule.edges(data=True):
        order = orders[str(data['bonds'])]
        rdkG.AddBond(int(u), int(v), order)

    rdkG = rdkG.GetMol()

    if sanitize:
        Chem.SanitizeMol(rdkG)

    return rdkG


def rdkit_to_gratoms(rdkG, name, valence, confid=0):
    block = Chem.MolToMolBlock(rdkG, confId=confid)

    positions = np.empty((rdkG.GetNumAtoms(), 3))
    symbols = []
    for i, atom in enumerate(block.split('\n')[4:rdkG.GetNumAtoms() + 4]):
        data = atom.split()
        positions[i] = np.array(data[:3], dtype=float)
        symbols += [data[3]]

    gratoms = Gratoms(symbols, positions)
    gratoms.graph.name = name
    nx.set_node_attributes(gratoms.graph, name='valence', values=valence)

    return gratoms


def plot_molecule(molecule, file_name=None):
    """Plot a molecule using RDKit."""
    rdkG = get_graph(molecule)
    rdkG = Chem.RemoveHs(rdkG)
    MolToFile(rdkG, file_name, size=(200, 200))


def get_smiles(molecule):
    """Return SMILES representation of a molecule as str."""
    rdkG = get_graph(molecule)

    return Chem.MolToSmiles(rdkG)


def get_uff_coordinates(gratoms, steps=10):
    rdkG = get_graph(gratoms.graph)
    Chem.EmbedMolecule(rdkG, Chem.ETKDG())

    lec = 0
    if steps:
        cids = Chem.EmbedMultipleConfs(rdkG, numConfs=steps)
        Chem.UFFOptimizeMoleculeConfs(rdkG)

        energies = []
        for cid in cids:
            ffe = Chem.UFFGetMoleculeForceField(rdkG, confId=cid).CalcEnergy()
            energies += [ffe]
        energies = np.array(energies)

        lec = int(np.argmin(energies))

    valence = nx.get_node_attributes(gratoms.graph, 'valence')
    name = gratoms.graph.name
    gratoms = rdkit_to_gratoms(rdkG, name, valence, confid=lec)

    return gratoms
