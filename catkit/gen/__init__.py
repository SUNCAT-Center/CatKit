# Copyright 2017
# (see accompanying license files for details).
"""Catalysis Generator."""

from collections.abc import MutableMapping
import numpy as np
import ase

radicals = np.ones(92)
radicals[[6, 7, 8, 9, 15, 16]] = [4, 3, 2, 1, 3, 2]


class Defaults(MutableMapping, dict):
    """No frills default dictionary class."""

    def __init__(self):
        self.update({
            'radii': ase.data.covalent_radii.copy(),
            'radicals': radicals,
            'orthogonal': False
        })

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

defaults = Defaults()

from . import symmetry
from . import adsorption
from . import surface
from . import molecules

__all__ = ['defaults', 'symmetry', 'adsorption', 'surface', 'molecules']
