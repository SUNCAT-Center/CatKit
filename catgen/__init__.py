# Copyright 2017
# (see accompanying license files for details).
"""Catalysis Generator."""

__version__ = '0.4.0'

from collections import MutableMapping
from ase.data import covalent_radii


class Defaults(MutableMapping, dict):
    """No frills default dictionary class."""

    def __init__(self):
        self.update({
            'covalent_radii': covalent_radii
        })

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)

    def __getitem__(self, key):
        return dict.__getitem__(self, key)


defaults = Defaults()

# Add a small skin to the default setting.
defaults['covalent_radii'] += 0.15
