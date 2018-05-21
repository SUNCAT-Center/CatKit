# Copyright 2017
# (see accompanying license files for details).
"""Catalysis Kit."""

from .gratoms import Gratoms
import matplotlib
matplotlib.use('Agg')

__all__ = ['Gratoms']
__version__ = '0.4.4'
