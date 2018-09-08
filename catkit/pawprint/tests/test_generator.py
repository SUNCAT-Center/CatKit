from catkit.pawprint import Fingerprinter
from catkit.build import surface
from ase.build import bulk
from ase.build import fcc111
from ase.build import add_adsorbate
import unittest
import numpy as np


class TestGenerator(unittest.TestCase):
    """Test features of catkit.pawprint.generator"""

    def test_nonlocal_fingerprinting(self):
        """Test the non-local fingerprint example"""
        atoms = bulk('Pd', cubic=True)
        atoms[3].symbol = 'Pt'

        slab = surface('Al', size=(2, 2, 3), a=3.8, vacuum=10)

        images = [atoms, slab]

        parameters = [
            'atomic_number',
            'atomic_radius',
            'atomic_volume',
        ]

        operations = [
            'periodic_convolution',
            ['periodic_convolution', {'d': 1}]
        ]

        fp = Fingerprinter(images)
        fingerprints = fp.get_fp(parameters, operations)

        truth = np.array([
            [12432.0, 7.562800000000001, 320.44, 136896.0, 90.7488, 3844.8],
            [2028.0, 24.53879999999999, 1200.0, 20280.0, 245.388, 12000.0]])

        np.testing.assert_allclose(fingerprints, truth)

    def test_local_fingerprinting(self):
        """Test the local fingerprint example"""
        atoms = fcc111('Pd', size=(2, 2, 3), vacuum=10.0)
        add_adsorbate(atoms, 'C', 1, 'fcc')

        # -1 is the tag convention to identify bonded atoms
        tags = atoms.get_tags()
        tags[-1] = -1
        atoms.set_tags(tags)

        parameters = [
            'atomic_number',
            'atomic_radius',
            'boiling_point',
            'covalent_radius_cordero',
            'evaporation_heat',
            'fusion_heat',
            'group_id',
            'period',
        ]

        operations = [
            'bonding_convolution',
            'layered_sum',
            'local_ads_metal_fp'
        ]

        fp = Fingerprinter(atoms)
        fingerprints = fp.get_fp(parameters, operations)

        truth = np.array([276.0, 1.2467000000000001, 17406300.0,
                          1.0147, np.nan, np.nan, 140.0, 10.0,
                          6.0000e+00, 1.8400e+02, 1.8400e+02, 1.8400e+02,
                          9.1000e-01, 5.4800e+00, 5.4800e+00, 5.4800e+00,
                          5.1000e+03, 1.3652e+04, 1.3652e+04, 1.3652e+04,
                          7.3000e-01, 5.5600e+00, 5.5600e+00, 5.5600e+00,
                              np.nan, 1.4896e+03, 1.4896e+03, 1.4896e+03,
                              np.nan, 6.8960e+01, 6.8960e+01, 6.8960e+01,
                          1.4000e+01, 4.0000e+01, 4.0000e+01, 4.0000e+01,
                          2.0000e+00, 2.0000e+01, 2.0000e+01, 2.0000e+01,
                          4.0000e+01, 4.6000e-01,-1.6870e+03, 6.6000e-01,
                              np.nan,     np.nan,-4.0000e+00, 3.0000e+00])

        np.testing.assert_allclose(fingerprints[0], truth)

    def test_custom_fingerprinting(self):
        """Test the custom fingerprint example"""
        atoms = fcc111('Pd', size=(2, 2, 3), vacuum=10.0)
        add_adsorbate(atoms, 'C', 1, 'fcc')

        # -1 is the tag convention to identify bonded atoms
        tags = atoms.get_tags()
        tags[-1] = -1
        atoms.set_tags(tags)

        parameters = [
            'atomic_number',
            'dband_center_slab',
            'dband_width_slab',
            'dband_skewness_slab',
            'dband_kurtosis_slab'
        ]

        def example_operation(
                atoms,
                atoms_parameters,
                connectivity):
            # This is a CatKit convention
            bond_index = np.where(atoms.get_tags() == -1)[0]

            return np.sum(connectivity[bond_index], axis=1)

        operations = [
            'bonding_convolution',
            example_operation
        ]

        fp = Fingerprinter(atoms)
        fingerprints = fp.get_fp(parameters, operations)

        truth = np.array([276.0, -1.570340288877117, 6.51684717398884,
                          -81.36785228863994, 3651.4867376228817, 3.0])

        np.testing.assert_allclose(fingerprints[0], truth)


if __name__ == '__main__':
    unittest.main()
