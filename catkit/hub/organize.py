#!/usr/bin/env python

# builtin imports
import pickle
import re
import pprint

# A lot of functions from os.path
# in python 2 moved to os. and changed
# their signature. Pathlib can be
# installed under python2.X with
# pip install pathlib2 and is in
# standard library in Python 3,
# hence we use it as a compatiblity
# library
try:
    from pathlib import Path
    Path().expanduser()
except (ImportError, AttributeError):
    from pathlib2 import Path


# other library imports
import Levenshtein
import ase.atoms
import ase.utils
import ase.io
import numpy as np
from ase.symbols import string2symbols


# local imports
from .ase_tools import gas_phase_references, get_chemical_formula, \
    symbols, collect_structures
import catkit.hub.ase_tools

np.set_printoptions(threshold=500, linewidth=1800, edgeitems=80)

PUBLICATION_TEMPLATE = """{
    "title": "Test",
    "authors": ["Doe"],
    "journal": "",
    "volume": "",
    "number": "",
    "pages": "",
    "year": "2018",
    "publisher": "",
    "doi": "",
    "tags": []
}"""


def fuzzy_match(structures, options):
    # filter out cell with ill-defined unit cells
    structures = [structure for structure in structures
                  if structure.number_of_lattice_vectors == 3
                  ]
    # sort by density
    structures = sorted(structures,
                        key=lambda x: len(x) / x.get_volume()
                        )

    # group in to bulk, surface, or bulk
    molecules, surfaces, bulks = [], [], []
    gas_phase_candidates = []
    reference_energy = {}
    collected_energies = {}
    key_count = {}
    collected_structures = {}
    if options.verbose:
        print("Group By Densities")
        print("===================")
    for structure in structures:
        if options.include_pattern:
            if not re.search(
                    options.include_pattern, structure.info['filename']):
                continue
            if options.exclude_pattern \
                    and re.search(
                        options.exclude_pattern,
                        structure.info['filename']):
                continue

        elif options.exclude_pattern:
            if re.search(
                    options.exclude_pattern,
                    structure.info['filename']):
                continue

        # add more info from filename
        facet_match = re.search(
            '(?<=[^0-9])?[0-9]{3,3}(?=[^0-9])?', structure.info['filename'])
        if facet_match:
            structure.info['facet'] = facet_match.group()
        else:
            structure.info['facet'] = options.facet_name or 'facet'

        density = len(structure) / structure.get_volume()
        if options.verbose:
            print("  {density:10.3f} {filename}".format(
                density=density,
                filename=structure.info['filename'],
            ))
        if density < options.max_density_gas:
            structure.info['state'] = 'molecule'
            molecules.append(structure)
            collected_structures \
                .setdefault(options.dft_code
                            or structure.info['filetype'], {}) \
                .setdefault(options.xc_functional or 'XC_FUNCTIONAL', {}) \
                .setdefault('gas', {}) \
                .setdefault(get_chemical_formula(structure), structure)

            formula = get_chemical_formula(structure)
            if formula not in options.exclude_reference.split(','):
                gas_phase_candidates.append(get_chemical_formula(structure))
                reference_energy[formula] = structure.get_potential_energy()
                if options.verbose:
                    print("           GAS", formula,
                          structure.info['filename'])

        elif density < options.max_density_slab:
            structure.info['state'] = 'surface'
            formula = get_chemical_formula(structure)
            surfaces.append(structure)
            if options.verbose:
                print("           SURFACE", formula,
                      structure.info['filename'])
        else:
            structure.info['state'] = 'bulk'
            bulks.append(structure)

            if options.verbose:
                print("           BULK", formula, structure.info['filename'])

    # sort surfaces by volume to get difference facets
    surfaces = sorted(surfaces,
                      key=lambda x: x.get_volume()
                      )

    # Get minimal set of gas phase candidates
    gas_phase_candidates = list(
        sorted(
            set(gas_phase_candidates)
        )
    )
    if options.verbose:
        print("\n\nCANDIDATES {gas_phase_candidates}".format(
            gas_phase_candidates=gas_phase_candidates,
        ))

    volume_groups = {}
    tolerance = 1e-5
    if options.verbose:
        print("\n\nGROUP BY VOLUME\n\n")
    for surface in sorted(surfaces,
                          key=lambda x: x.get_volume(),):
        formula = symbols(surface)

        for volume in volume_groups:
            if abs(volume - surface.get_volume()) < tolerance:
                volume_groups[volume].append(surface)
                break
        else:
            volume_groups[surface.get_volume()] = [surface]
            if options.verbose:
                print(("\n=============== NEW VOLUME"
                       " {volume} ============"
                       ).format(volume=surface.get_volume()))
        if options.verbose:
            print(get_chemical_formula(surface))

    for volume in volume_groups:
        if options.verbose:
            print("\nInspect volume {volume}\n".format(
                volume=volume,
            ))
        surfaces = volume_groups[volume]
        N = len(surfaces)
        if N > 1:
            distance = np.zeros((N, N), dtype=int)
            distance[:] = 99
            for i, structure1 in enumerate(surfaces):
                formula1 = symbols(structure1)
                for j, structure2 in enumerate(surfaces):
                    if j >= i:
                        formula2 = symbols(structure2)
                        distance[i, j] = Levenshtein.distance(
                            formula1, formula2)

        for i, surf1 in enumerate(surfaces):
            for j, surf2 in enumerate(surfaces):
                f1 = symbols(surf1)
                formula1 = get_chemical_formula(surf1)
                f2 = symbols(surf2)
                formula2 = get_chemical_formula(surf2)
                if Levenshtein.distance(f1, f2) in range(1, 10):
                    additions = ''
                    subtractions = ''
                    equal = ''
                    opcodes = Levenshtein.opcodes(f1, f2)
                    for tag, i1, i2, j1, j2 in opcodes:
                        if tag == 'insert':
                            additions += f2[j1:j2]
                        elif tag == 'replace':
                            additions += f2[j1:j2]
                        elif tag == 'delete':
                            subtractions += f2[j1:j2]
                        elif tag == 'equal':
                            equal += f1[i1:i2]

                    if options.verbose:
                        print("    ADDITIONS " + str(additions))
                        print("    SUBTRACTIONS " + str(subtractions))

                    try:
                        subtractions = ''.join(
                            sorted(
                                string2symbols(
                                    subtractions)))
                    except Exception as e:
                        if options.verbose:
                            print("Warning: trouble parsing {subtractions}:{e}"
                                  .format(subtractions=subtractions, e=e))
                    try:
                        additions = ''.join(
                            sorted(
                                    string2symbols(
                                        additions)))
                    except Exception as e:
                        if options.verbose:
                            print("Warning: trouble parsing {additions}, {e}"
                                  .format(additions=additions, e=e))

                    equal_formula = get_chemical_formula(
                        ase.atoms.Atoms(equal))

                    # check if we have some additions of subtractions
                    # and either one (or both) are in user specifid
                    # adsorbates
                    if (additions or subtractions) \
                            and (not additions or
                                 additions in options.adsorbates) \
                            and (not subtractions or
                                 subtractions in options.adsorbates):

                        dE = surf2.get_potential_energy() \
                            - surf1.get_potential_energy()
                        difference_symbols = \
                            gas_phase_references.molecules2symbols(
                                [additions, subtractions],
                            )

                        adsorbates = []
                        if additions:
                            adsorbates.append(additions)
                        if subtractions:
                            adsorbates.append(subtractions)

                        if options.verbose:
                            print("    ADDITIONS " + str(additions))
                            print("    SUBTRACTIONS " + str(subtractions))
                            print("    ADSORBATES " + str(adsorbates))

                        # TODO: len(gas_phase_candidates) >= symbols
                        if len(gas_phase_candidates)  \
                           >= len(difference_symbols):
                            references = \
                                gas_phase_references \
                                .construct_reference_system(
                                    difference_symbols,
                                    gas_phase_candidates,
                                    options,
                                )
                            if options.verbose:
                                print(" REFERENCES " + str(references))

                            stoichiometry_factors =  \
                                gas_phase_references.get_stoichiometry_factors(
                                    adsorbates, references,
                                )
                        else:
                            adsorbates = map(lambda x: ase.utils.formula_hill(
                                            catkit.hub.ase_tools.get_numbers_from_formula(x))
                                    , adsorbates)
                            stoichiometry_factors = {}
                            if options.verbose:
                                print(" ADSORBATES " + str(adsorbates))
                                print(" GP_CANDIDATES " + str(gas_phase_candidates))
                            for adsorbate in adsorbates:
                                if adsorbate in gas_phase_candidates:
                                    stoichiometry_factors \
                                        .setdefault(adsorbate, {}) \
                                        .setdefault(adsorbate, 1)
                                else:
                                    raise UserWarning((
                                        "Could not construct stoichiometry"
                                        " factors for {adsorbate}\n"
                                        "from {candidates}."
                                        "Please add more gas phase molecules"
                                        " to your folder.\n"
                                    ).format(
                                        adsorbate=adsorbate,
                                        candidates=gas_phase_candidates,
                                    ))

                        if options.verbose:
                            print("STOICHIOMETRY FACTORS "
                                  + str(stoichiometry_factors))
                        if options.verbose:
                            print("COLLECTED ENERGIES")
                            print(collected_energies)
                            print("    STOICH FACTORS " +
                                  str(stoichiometry_factors) + "\n\n")

                        matching_keys = [
                            _key for _key in collected_energies
                            if _key.startswith(key)
                        ]

                        adsorbate = get_chemical_formula(
                            ase.atoms.Atoms(additions))

                        key = ("{equal_formula}"
                               "({surface_facet})"
                               "+{adsorbate}"
                               ).format(
                            formula=formula,
                            equal_formula=equal_formula,
                            surface_facet=surf1.info['facet'],
                            adsorbate=adsorbate,
                        )

                        formula = '*@site' + str(key_count.get(key, 0)) + ' ->'

                        if additions:
                            formula += ' ' + \
                                get_chemical_formula(
                                    ase.atoms.Atoms(additions)) \
                                + '*@site' + str(key_count.get(key, 0))
                        if subtractions:
                            formula += ' ' + \
                                get_chemical_formula(
                                    ase.atoms.Atoms(subtractions)) \
                                + '*@site' + str(key_count.get(key, 0))

                        gas_phase_corrections = {}

                        for adsorbate in adsorbates:
                            stoich_factors = stoichiometry_factors[adsorbate]
                            for ref in stoich_factors:
                                dE -= stoich_factors[ref] * \
                                    reference_energy[ref]
                                gas_phase_corrections[ref] = \
                                    gas_phase_corrections.get(
                                        ref, 0) - stoich_factors[ref]

                        for molecule, factor in gas_phase_corrections.items():
                            if factor != 0:
                                sign = ' + ' if factor < 0 else ' +- '
                                if abs(factor - int(factor)) < 1e-3:
                                    factor = str(abs(int(factor)))
                                    if factor == '1':
                                        factor = ''
                                else:
                                    factor = '{:.2f}'.format(abs(factor))

                                fleft, fright = formula.split(' -> ')
                                formula = fleft + sign + factor + \
                                    molecule + '(g)' + ' -> ' + fright

                        if abs(dE) < options.max_energy:
                            energy = dE
                            if options.verbose:
                                print("KEY {key}".format(**locals()))
                            equation = (" {formula:30s}"
                                        ).format(
                                formula=formula,
                                equal_formula=equal_formula,
                                surface_facet=surf1.info['facet']
                            ) \
                                .replace(' ', '') \
                                .replace('+', '_') \
                                .replace('->', '__') \
                                .replace('*', 'star') \
                                .replace('(g)', 'gas')
                            # We keep the empty structure whether or not
                            # we keep all structures
                            collected_structures \
                                .setdefault(
                                    options.dft_code
                                    or structure.info['filetype'],
                                    {}) \
                                .setdefault(options.xc_functional, {}) \
                                .setdefault(equal_formula + '_' + (
                                    options.structure
                                    or 'structure'
                                    ), {}) \
                                .setdefault(
                                        options.facet_name
                                        if options.facet_name != 'facet'
                                        else surf1.info['facet']
                                        , {}) \
                                .setdefault('empty_slab', surf1)

                            collected_energies[key] = energy
                            key_count[key] = key_count.get(key, 0) + 1
                            if options.verbose:
                                print(key)
                                print(collected_energies)
                                print(key in collected_energies)
                            if not options.keep_all_energies:
                                if energy > collected_energies.get(
                                        key, float("inf")):
                                    continue

                            # persist adsorbate slab structures
                            ####################################
                            collected_energies[key] = energy
                            collected_structures .setdefault(
                                options.dft_code
                                or structure.info['filetype'],
                                {}) .setdefault(
                                options.xc_functional,
                                {}) .setdefault(
                                equal_formula + '_' + (
                                    options.structure
                                    or 'structure'), {}
                                ).setdefault(
                                options.facet_name
                                if options.facet_name != 'facet'
                                else surf1.info['facet'],
                                {}) .setdefault(
                                equation,
                                {})[adsorbate] = surf2

    print("\n\nCollected Reaction Energies Data")
    print("====================================")
    if options.verbose:
        pprint.pprint(collected_structures)
    print("\n\nCollected Reaction Energies")
    print("===========================")
    if len(collected_energies) == 0:
        print("Warning: no energies collected. Some ways to fix this:")
        print("  * raise the allowed maximum reaction energy (default: 10 eV)")
        print("     --max-energy 50")
        print("  * make sure you have gas phase molecules in the directory")
        print("  * raise the maximum density for gas-phase molecule")
        print("     --max-density-gas 0.004")
        print("  * raise the maximum density for slab structures")
        print("    --max-density-slab 0.03 ")

    for key, energy in collected_energies.items():
        print("{key:40s}: {energy:.3f} eV".format(
            key=key,
            energy=energy,
        ))

    return collected_structures


def create_folders(options, structures, root='',
                   publication_template=PUBLICATION_TEMPLATE):
    out_format = 'json'

    for key in structures:
        if isinstance(structures[key], dict):
            d = Path(root).joinpath(key)
            Path(d).mkdir(parents=True, exist_ok=True)
            if Path(root).parent.as_posix() == '.':
                # Have to explicitly convert Path to str
                # to work under python 3.4
                with open(str(
                        Path(root).joinpath('publication.txt')),
                        'w') as outfile:
                    outfile.write(publication_template)
            create_folders(options, structures[key], root=d)
        else:
            ase.io.write(
                str(Path(root).joinpath(key + '.' + out_format)),
                structures[key],
                format=out_format,
            )


def main(options):
    pickle_file = options.foldername.strip().rstrip(
        '/').strip('.').rstrip('/') + '.cache.pckl'

    if Path(pickle_file).exists() \
            and Path(pickle_file).stat().st_size \
            and options.use_cache:
        with open(pickle_file, 'rb') as infile:
            structures = pickle.load(infile)
    else:
        structures = collect_structures(options.foldername, options.verbose,
                                        level='**/*')
        if options.gas_dir:
            structures.extend(
                collect_structures(
                    options.gas_dir,
                    options.verbose,
                    level='**/*')
                    )
        if options.use_cache:
            with open(pickle_file, 'wb') as outfile:
                pickle.dump(structures, outfile)

    if hasattr(catkit.hub.ase_tools, 'PUBLICATION_TEMPLATE') \
            and catkit.hub.ase_tools.PUBLICATION_TEMPLATE:
        publication_template = catkit.hub.ase_tools.PUBLICATION_TEMPLATE
    else:
        publication_template = PUBLICATION_TEMPLATE

    structures = fuzzy_match(structures, options)
    create_folders(options, structures,
                   root=options.foldername.strip('/') + '.organized',
                   publication_template=publication_template,
                   )
