# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 21:01:59 2018

@author: mhangaard
"""
import numpy as np
from ase.build import bulk
from catkit.surface import SlabGenerator
from ase.build import add_adsorbate
from ase.data import covalent_radii
from atoml.fingerprint.database_adsorbate_api import (get_radius,
                                                      attach_adsorbate_info)
from atoml.fingerprint.setup import FeatureGenerator
from atoml.preprocess.clean_data import clean_infinite


# Define an adsorbate.
ads = 'C'

# Create a bulk.
atoms = bulk('Pd', 'fcc', a=4, cubic=True)
atoms[3].symbol = 'Cu'

# Create some surfaces.
Gen = SlabGenerator(
    atoms,
    miller_index=[1, 1, 1],
    layers=3,
    fixed=1,
    vacuum=8,
)
terminations = Gen.get_unique_terminations()
# Loop over unique surface terminations.
images = []
for i, t in enumerate(terminations):
    slab = Gen.get_slab(iterm=i)
    slab.center(axis=2, vacuum=5)
    # Add an adsorbate to each slab.
    radii = np.average([get_radius(a.number) for a in slab])
    height = radii + covalent_radii[6]
    add_adsorbate(slab, ads, height, position=slab.get_positions()[-1, :2])
    # Some help is needed for AtoML to positively identify the adsorbate.
    # A string in Hill notation is enough, if the adsorbate consists
    #   of different elements from the slab.
    slab.info['species'] = ads
    # The number of layers also makes some of the fingerprinting easier.
    slab.info['layers'] = Gen.layers
    # Append them to a list.
    images.append(slab)

# Some information is expected to be attached to atoms objects.
#   There are various ways of doing this, but the easiest is to call
AtoML_atoms = attach_adsorbate_info(images)
# This is where checks should be made

# Instantiate the fingerprint generator for adsorbate structures.
fingerprinter = FeatureGenerator()

# All user methods under the fingerprinter accepts an atoms object and
#   returns a vector.
functions = [
    fingerprinter.ads_nbonds,
    fingerprinter.primary_addatom,
    fingerprinter.primary_adds_nn,
    fingerprinter.Z_add,
    fingerprinter.ads_av,
    fingerprinter.primary_surf_nn,
    fingerprinter.primary_surfatom,
    fingerprinter.primary_surfatom_sum,
    ]

# This list is passed on to the following setup functions,
#    along with a list of atoms.

# Get and print the names of features.
features_labels = fingerprinter.return_names(functions)
for l in range(len(features_labels)):
    print(l, features_labels[l])

# Get a matrix containing the fingerprints.
unlabeled_data_matrix = fingerprinter.return_vec(AtoML_atoms, functions)
print(np.shape(unlabeled_data_matrix), 'data matrix created.')

# Cleanup in case some of the functions are returning NaNs or Infs
clean_data_matrix = clean_infinite(unlabeled_data_matrix)['train']

# Ready for Machine learning.
print(np.shape(clean_data_matrix), 'data matrix returned.')
