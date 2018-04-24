from catkit import Gratoms
from ase.data import chemical_symbols
import numpy as np
import networkx as nx
from rdkit.Chem.Draw import MolToFile
from rdkit.Chem import AllChem as Chem


def get_graph(molecule, sanitize=True):

    rdkG = Chem.rdchem.EditableMol(Chem.rdchem.Mol())

    for j, data in molecule.nodes(data=True):
        rdAtom = Chem.rdchem.Atom(chemical_symbols[data['number']])
        if data.get('valence'):
            rdAtom.SetNumRadicalElectrons(int(data.get('valence')))
        else:
            rdAtom.SetNumRadicalElectrons(int(0))
        rdkG.AddAtom(rdAtom)

    rdBonds = Chem.rdchem.BondType
    orders = {'1': rdBonds.SINGLE, '2': rdBonds.DOUBLE, '3': rdBonds.TRIPLE}

    for u, v, data in molecule.edges(data=True):
        if data.get('bonds'):
            order = orders[str(data['bonds'])]
        else:
            order = rdBonds.SINGLE
        rdkG.AddBond(int(u), int(v), order)

    rdkG = rdkG.GetMol()

    if sanitize:
        Chem.SanitizeMol(rdkG)

    return rdkG


def rdkit_to_gratoms(rdkG, name=None, confid=-1):
    """TODO: conserve 3D positions if present."""
    block = Chem.MolToMolBlock(rdkG, confId=confid)

    positions = np.empty((rdkG.GetNumAtoms(), 3))
    n = rdkG.GetNumAtoms()

    symbols, edges, valence = [], [], {}
    for i, atom in enumerate(block.split('\n')[4:n + 4]):
        data = atom.split()
        positions[i] = np.array(data[:3], dtype=float)
        symbols += [data[3]]
        valence.update({i: int(data[9])})

    for i, bond in enumerate(block.split('\n')[n + 4:]):
        data = bond.split()
        if data[0] == 'M':
            break
        data = np.array(data, dtype=int)

        edges += [(data[0] - 1, data[1] - 1, {'bonds': data[2]})]

    gratoms = Gratoms(symbols, positions)
    gratoms.graph.name = name
    nx.set_node_attributes(gratoms.graph, values=valence, name='valence')
    gratoms.graph.add_edges_from(edges)

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

    name = gratoms.graph.name
    gratoms = rdkit_to_gratoms(rdkG, name, confid=lec)

    return gratoms
