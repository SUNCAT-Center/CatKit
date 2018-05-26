#!/usr/bin/env python

import subprocess
import sys
import os
import collections

import catkit.hub.organize

print(os.path.abspath(__file__))
path = os.path.abspath(os.path.join(os.path.dirname(__file__)))


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

# def test_debugging():
    #raise UserWarning(sys.path)


def test_file_organization():

    subprocess.call(
        ('python {path}/make_test_slabs.py'.format(path=path)).split())
    subprocess.call(
        ('cathub organize {path}/unorganized --adsorbates O,H2'.format(path=path)).split())


def test_file_organization_module():
    options = Struct(**{
        'adsorbates': 'O,H2',
        'foldername': '{path}/unorganized'.format(path=path),
        'verbose': True,
        'include_pattern': '.',
        'exclude_pattern': '%%$^#$',
        'facet_name': '111',
        'max_density_gas': 0.002,
        'max_density_slab': 0.08,
        'exclude_reference': '',
        'max_energy': 10,
        'keep_all_energies': True,

    })

    catkit.hub.organize.main(options)


if __name__ == '__main__':
    test_file_organization()
    test_file_organization_module()
