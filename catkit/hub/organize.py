#!/usr/bin/env python

# builtin imports
import pickle
import re

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
import ase.io
import numpy as np


# local imports
from .ase_tools import gas_phase_references

np.set_printoptions(threshold=500, linewidth=1800, edgeitems=80)

PUBLICATION_TEMPLATE = """{
    "title": "",
    "authors": [""],
    "journal": "",
    "volume": "",
    "number": "",
     "pages": "",
     "year": "",
     "publisher": "",
     "doi": "",
     "tags": []
}"""


def get_chemical_formula(atoms):
    """
    Compatibility function, return mode=metal, when
    available, mode=hill, when not (ASE <= 3.13)
    """
    try:
        return atoms.get_chemical_formula(mode='metal')
    except ValueError:
        return atoms.get_chemical_formula(mode='hill')


def symbols(atoms):
    formula = get_chemical_formula(atoms)
    symbols = ase.atoms.string2symbols(formula)
    return ''.join(symbols)


def collect_structures(foldername, options):
    structures = []
    if options.verbose:
        print(foldername)
    for i, filename in enumerate(Path(foldername).glob('**/*')):
        posix_filename = str(filename.as_posix())
        if options.verbose:
            print(i, posix_filename)
        if posix_filename.endswith('publication.txt'):
            with open(posix_filename) as infile:

                global PUBLICATION_TEMPLATE
                PUBLICATION_TEMPLATE = infile.read()
        elif Path(posix_filename).is_file():
            try:
                filetype = ase.io.formats.filetype(posix_filename)
            except Exception as e:
                continue
            if filetype:
                try:
                    structure = ase.io.read(posix_filename)
                    structure.info['filename'] = posix_filename
                    structure.info['filetype'] = ase.io.formats.filetype(
                        posix_filename)
                    try:
                        structure.get_potential_energy()
                        # ensure that the structure has an energy
                        structures.append(structure)
                    except RuntimeError:
                        print("Did not add {posix_filename} since it has no energy"
                              .format(
                                  posix_filename=posix_filename,
                              ))
                    print(structure)
                except StopIteration:
                    print("Warning: StopIteration {posix_filename} hit."
                          .format(
                              posix_filename=posix_filename,
                          ))
                except IndexError:
                    print("Warning: File {posix_filename} looks incomplete"
                          .format(
                              posix_filename=posix_filename,
                          ))
                except OSError as e:
                    print("Error with {posix_filename}: {e}".format(
                        posix_filename=posix_filename,
                        e=e,
                    ))
                except AssertionError as e:
                    print("Hit an assertion error with {posix_filename}: {e}".format(
                        posix_filename=posix_filename,
                        e=e,
                    ))
                except ValueError as e:
                    print("Trouble reading {posix_filename}: {e}".format(
                        posix_filename=posix_filename,
                        e=e,
                    ))

    return structures


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
                .setdefault(structure.info['filetype'], {}) \
                .setdefault('XC_FUNCTIONAL', {}) \
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
                        elif tag == 'delete':
                            subtractions += f2[j1:j2]
                        elif tag == 'equal':
                            equal += f1[i1:i2]
                    subtractions = ''.join(
                        sorted(
                            ase.atoms.string2symbols(
                                subtractions)))
                    additions = ''.join(
                        sorted(
                                ase.atoms.string2symbols(
                                    additions)))

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
                        references = \
                            gas_phase_references.construct_reference_system(
                                difference_symbols, gas_phase_candidates,
                            )
                        atomic_references = {}
                        for i, j in references:
                            atomic_references[i] = j

                        adsorbates = []
                        if additions:
                            adsorbates.append(additions)
                        if subtractions:
                            adsorbates.append(subtractions)
                        if options.verbose:
                            print("    ADDITIONS " + str(additions))
                            print("    SUBTRACTIONS " + str(subtractions))
                            print("    ADSORBATES " + str(adsorbates))
                            print("    REFERENCES " + str(references))
                        stoichiometry_factors =  \
                            gas_phase_references.get_stoichiometry_factors(
                                adsorbates, references,
                            )
                        if options.verbose:
                            print("    STOICH FACTORS " +
                                  str(stoichiometry_factors) + "\n\n")

                        formula = '* ->'
                        adsorbate = get_chemical_formula(
                            ase.atoms.Atoms(additions))
                        if additions:
                            formula += ' ' + \
                                get_chemical_formula(
                                    ase.atoms.Atoms(additions)) + '*'
                        if subtractions:
                            formula += ' ' + \
                                get_chemical_formula(
                                    ase.atoms.Atoms(subtractions)) + '*'

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

                            key = ("{equal_formula:16s}"
                                   " ({surface_facet})"
                                   " {formula:30s}"
                                   ).format(
                                formula=formula,
                                equal_formula=equal_formula,
                                surface_facet=surf1.info['facet']
                            )
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
                                .setdefault(structure.info['filetype'], {}) \
                                .setdefault('XC_FUNCTIONAL', {}) \
                                .setdefault(equal_formula, {}) \
                                .setdefault(surf1.info['facet'], {}) \
                                .setdefault('empty_slab', surf1)
                            if not options.keep_all_energies:
                                if energy < collected_energies.get(
                                        key, float("inf")):
                                    collected_energies[key] = energy
                                    collected_structures .setdefault(
                                        structure.info['filetype'],
                                        {}) .setdefault(
                                        'XC_FUNCTIONAL',
                                        {}) .setdefault(
                                        equal_formula,
                                        {}) .setdefault(
                                        surf1.info['facet'],
                                        {}) .setdefault(
                                        equation,
                                        {})[adsorbate] = surf2
                            else:
                                # amend number if colliding keys are found
                                matching_keys = [
                                    _key for _key in collected_energies
                                    if _key.startswith(key)
                                ]
                                key += '_' + str(len(matching_keys))
                                collected_energies[key] = energy
                                collected_structures .setdefault(
                                    structure.info['filetype'], {}) \
                                    .setdefault('XC_FUNCTIONAL', {}) \
                                    .setdefault(equal_formula, {}) \
                                    .setdefault(surf1.info['facet'], {}) \
                                    .setdefault(
                                    equation, {}) .setdefault(
                                    adsorbate + '@site_' + str(
                                        len(matching_keys)), surf2)

    print("\n\nCollected Reaction Energies")
    print("===========================")
    for key, energy in collected_energies.items():
        print("{key:40s}: {energy:.3f} eV".format(
            key=key,
            energy=energy,
        ))

    return collected_structures


def create_folders(options, structures, root=''):
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
                    outfile.write(PUBLICATION_TEMPLATE)
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

    if Path(pickle_file).exists() and Path(pickle_file).stat().st_size:
        with open(pickle_file, 'rb') as infile:
            structures = pickle.load(infile)
    else:
        structures = collect_structures(options.foldername, options)
        with open(pickle_file, 'wb') as outfile:
            pickle.dump(structures, outfile)

    structures = fuzzy_match(structures, options)
    create_folders(options, structures,
                   root=options.foldername.strip('/') + '.organized')
