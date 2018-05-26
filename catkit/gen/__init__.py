# Copyright 2017
# (see accompanying license files for details).
"""Catalysis Generator."""

from collections import MutableMapping
from ase.data import covalent_radii
import numpy as np

radicals = np.ones(92)
radicals[[6, 7, 8, 9, 15, 16]] = [4, 3, 2, 1, 3, 2]


class Defaults(MutableMapping, dict):
    """No frills default dictionary class."""

    def __init__(self):
        self.update({
            'covalent_radii': covalent_radii.copy(),
            'radicals': radicals,
        })

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)

    def __getitem__(self, key):
        return dict.__getitem__(self, key)


defaults = Defaults()
