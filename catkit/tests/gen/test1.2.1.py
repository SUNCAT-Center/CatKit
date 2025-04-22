from catkit.gen.adsorption import AdsorptionSites
from catkit.gen.surface import SlabGenerator
from ase.build import bulk

bulk = bulk('Pd', 'fcc', a=5, cubic=True)
bulk[3].symbol = 'Cu'

gen = SlabGenerator(
    bulk,
    miller_index=(1, 1, 1),
    layers=6,
    vacuum=4)

atoms = gen.get_slab()
atoms.set_surface_atoms([8, 9, 10, 11])

sites = AdsorptionSites(atoms)
sites.plot('./Pd3Cu-adsorption-sites.png')
