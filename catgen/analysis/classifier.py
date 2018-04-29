from .. import utils
import numpy as np


class Classifier():
    """Class for classification of various aspects of an an atomic
    unit cell.

    Currently, a tool for classification of adsorbates on surface
    environments and the active sites they rest on.
    """

    def __init__(self, atoms):
        """Return unique coordinate values of a given atoms object
        for a specified axis.

        Parameters:
        -----------
            atoms : atoms object
        """
        self.atoms = atoms
        self.ads_atoms = None
        self.slab_atoms = None

    def id_slab_atoms(
            self,
            classifier='trivial',
            tag=False,
            rtol=1e-3):
        """Return the indices of the slab atoms using select characterization
        techniques.

        Parameters:
        -----------
        classifier : str
            Classification technique to identify slab atoms.

            'trivial':
            Slab atoms assumed to have atomic number == 13 or >= 21.
        tag : bool
            Return adsorbate atoms with tags of 2.
        rtol : float
            Relative cutoff distance for tagging layers.

        Returns:
        --------
        slab_atoms : ndarray (n,)
            Index of slab atoms found.
        """
        atoms = self.atoms

        if classifier == 'trivial':
            slab_atoms = np.where((atoms.numbers == 13) |
                                  (atoms.numbers >= 21))[0]

        if tag:
            zpos = np.sort(atoms.positions[slab_atoms][:, -1])
            new_tags = np.zeros_like(zpos, dtype=int)
            tag = 1
            for i, z in enumerate(zpos):
                if new_tags[i] != 0:
                    continue
                layer = np.isclose(z, zpos, rtol=rtol)
                new_tags[layer] = tag
                tag += 1

            tags = self.atoms.get_tags()
            tags[slab_atoms] = new_tags[::-1]
            self.atoms.set_tags(tags)

        return slab_atoms

    def id_adsorbate_atoms(self, classifier='trivial', tag=False):
        """Identify adsorbed atoms in a given atoms object.

        Parameters:
        -----------
        classifier : str
            Classification technique to identify adsorbate atoms.

            'trivial':
            Adsorbate atoms assumed to have atomic number != 13 or < 21.
        tag : bool
           Return adsorbate atoms with tags of -2.

        Returns:
        --------
        ads_atoms : ndarray (n,)
            Index of adsorbate atoms found.
        """
        atoms = self.atoms

        if classifier == 'trivial':
            ads_atoms = np.where((atoms.numbers != 13) &
                                 (atoms.numbers < 21))[0]

        if tag:
            tags = self.atoms.get_tags()
            tags[ads_atoms] = -2
            self.atoms.set_tags(tags)

        self.ads_atoms = ads_atoms

        return ads_atoms
