#!/usr/bin/env python

import copy

import ase
import ase.build
import ase.lattice
import ase.io
import ase.lattice.surface
import ase.calculators.singlepoint


slab1 = ase.lattice.surface.fcc111('Pt', [2, 2, 4], vacuum=10)
slab1.set_calculator(
        ase.calculators.singlepoint.SinglePointCalculator(
            atoms=slab1,
            energy=500,
            )
        )
ase.io.write('empty_slab_111.traj', slab1)
ase.lattice.surface.add_adsorbate(slab1, ase.atoms.Atoms('O'), height=1.5)
slab1.set_calculator(
        ase.calculators.singlepoint.SinglePointCalculator(
            atoms=slab1,
            energy=510,
            )
        )
ase.io.write('empty_slab_111_ads.traj', slab1)

mol = ase.atoms.Atoms('O',
        cell=[
            [15., 0., 0.],
            [0., 15., 0.],
            [0, 0., 15.],
            ]
            )
mol.set_calculator(
        ase.calculators.singlepoint.SinglePointCalculator(
            energy=13,
            atoms=mol,
            )
        )
ase.io.write('ads.traj', mol )

mol = ase.atoms.Atoms('H2',
        cell=[
            [15., 0., 0.],
            [0., 15., 0.],
            [0, 0., 15.],
            ]
            )
mol.set_calculator(
        ase.calculators.singlepoint.SinglePointCalculator(
            energy=0.13,
            atoms=mol,
            )
        )
ase.io.write('ads2.traj', mol )
