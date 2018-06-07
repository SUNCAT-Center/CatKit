import ase.atoms
import ase.data
import numpy as np
import optparse
import pprint


def molecules2symbols(molecules, add_hydrogen=True):
    """Take a list of molecules and return just a list of atomic
    symbols, possibly adding hydrogen
    """
    symbols = sorted(
        list(set(
            ase.atoms.string2symbols(''.join(
                map(
                    lambda _x:
                    ''.join(ase.atoms.string2symbols(_x)), molecules)
            ))
        )),
        key=lambda _y: ase.data.atomic_numbers[_y])

    if add_hydrogen and 'H' not in symbols:
        symbols.insert(0, 'H')

    return symbols


def construct_reference_system(
    symbols,
    candidates=[
        'H2',
        'H2O',
        'NH3',
        'CH4',
        'CO',
        'SH2',
        'HCl',
        'N2',
        'O2']):
    """Take a list of symbols and construct gas phase
    references system, when possible avoiding O2.
    Candidates can be rearranged, where earlier candidates
    get higher preference than later candidates

    assume symbols sorted by atomic number
    """
    references = {}
    added_symbols = []
    for symbol in symbols:
        added_symbols.append(symbol)
        for candidate in candidates:
            symbols = ase.atoms.string2symbols(candidate)

            if set(added_symbols) == set(list(references.keys()) + symbols):
                references[symbol] = candidate
                break
        else:
            raise UserWarning(
                "No candidate satisfied {symbol}. Add more candidates".format(
                    **locals()))

    return sorted(references.items(),
                  key=lambda _x: ase.data.atomic_numbers[_x[0]]
                  )


def get_atomic_stoichiometry(references):
    """Given a list of references (tuples of (symbol, molecule))
    return stoichiometry matrix that connects atomic symbols
    to its molecular references.
    """
    n = len(references)
    stoichiometry = np.ndarray((n, n), float)
    stoichiometry[:] = 0.
    key_index = {}
    for i, (key, species) in enumerate(references):

        # in case species uses a suffix like _gas
        species = species.split('_')[0]

        key_index[key] = i
        composition = ase.atoms.string2symbols(species)
        for j, symbol in enumerate(composition):
            if symbol == key:
                stoichiometry[i, i] += 1
            else:
                stoichiometry[i, key_index[symbol]] -= 1
    for i in range(n):
        stoichiometry[i, i] = 1. / stoichiometry[i, i]

    for i in range(n):
        for j in range(n):
            if i != j:
                stoichiometry[i, j] *= stoichiometry[j, j]

    return stoichiometry.tolist()


def get_stoichiometry_factors(adsorbates, references):
    """Take a list of adsorabtes and a corresponding reference
    system and return a list of dictionaries encoding the
    stoichiometry factors converting between adsorbates and
    reference molecules.
    """
    stoichiometry = get_atomic_stoichiometry(references)
    stoichiometry_factors = {}
    for adsorbate in adsorbates:
        for symbol in ase.atoms.string2symbols(adsorbate):
            symbol_index = list(
                map(lambda _x: _x[0], references)).index(symbol)

            for (factor, (ref_symbol, ref_molecule)) in zip(
                    stoichiometry[symbol_index], references):
                stoichiometry_factors.setdefault(
                    adsorbate,
                    {})[ref_molecule] = stoichiometry_factors.setdefault(
                    adsorbate,
                    {}).get(
                    ref_molecule,
                    0) + factor

        nonzero_factors = {}
        for key, value in stoichiometry_factors[adsorbate].items():
            if not np.isclose(value, 0.):
                nonzero_factors[key] = value
        stoichiometry_factors[adsorbate] = nonzero_factors

    return stoichiometry_factors


if __name__ == '__main__':

    # store previously create test results
    test_results = [{'adsorbates': ['NH', 'CO', 'O', 'SH', 'OH', 'CH3'],
                     'references': [('H', 'H2'),
                                    ('C', 'CH4'),
                                    ('N', 'NH3'),
                                    ('O', 'H2O'),
                                    ('S', 'H2S')],
                     'stoichiometry': [[0.5, 0.0, 0.0, 0.0, 0.0],
                                       [-2.0, 1.0, 0.0, 0.0, 0.0],
                                       [-1.5, 0.0, 1.0, 0.0, 0.0],
                                       [-1.0, 0.0, 0.0, 1.0, 0.0],
                                       [-1.0, 0.0, 0.0, 0.0, 1.0]],
                     'stoichiometry_factors': {'CH3': {'CH4': 1.0, 'H2': -0.5},
                                               'CO': {'CH4': 1.0,
                                                      'H2': -3.0,
                                                      'H2O': 1.0},
                                               'NH': {'H2': -1.0,
                                                      'NH3': 1.0},
                                               'O': {'H2': -1.0,
                                                     'H2O': 1.0},
                                               'OH': {'H2': -0.5,
                                                      'H2O': 1.0},
                                               'SH': {'H2': -0.5,
                                                      'H2S': 1.0}},
                     'symbols': ['H', 'C', 'N', 'O', 'S']},
                    {'adsorbates': ['CO', 'O'],
                     'references': [('H', 'H2'), ('C', 'CH4'), ('O', 'H2O')],
                     'stoichiometry': [[0.5, 0.0, 0.0],
                                       [-2.0, 1.0, 0.0],
                                       [-1.0, 0.0, 1.0]],
                     'stoichiometry_factors': {'CO':
                                               {'CH4': 1.0,
                                                'H2': -3.0,
                                                'H2O': 1.0},
                                               'O': {'H2': -1.0, 'H2O': 1.0}},
                     'symbols': ['H', 'C', 'O']},
                    {'adsorbates': ['CO', 'N'],
                     'references': [
                         ('H', 'H2'),
                         ('C', 'CH4'),
                         ('N', 'NH3'),
                         ('O', 'H2O')],
                     'stoichiometry': [[0.5, 0.0, 0.0, 0.0],
                                       [-2.0, 1.0, 0.0, 0.0],
                                       [-1.5, 0.0, 1.0, 0.0],
                                       [-1.0, 0.0, 0.0, 1.0]],
                     'stoichiometry_factors': {
                         'CO': {'CH4': 1.0, 'H2': -3.0, 'H2O': 1.0},
                        'N': {'H2': -1.5, 'NH3': 1.0}},
                     'symbols': ['H', 'C', 'N', 'O']},
                    {'adsorbates': ['NO', 'O2', 'CO', 'CO2'],
                     'references': [('H', 'H2'),
                                    ('C', 'CH4'),
                                    ('N', 'NH3'),
                                    ('O', 'H2O')],
                     'stoichiometry': [[0.5, 0.0, 0.0, 0.0],
                                       [-2.0, 1.0, 0.0, 0.0],
                                       [-1.5, 0.0, 1.0, 0.0],
                                       [-1.0, 0.0, 0.0, 1.0]],
                     'stoichiometry_factors': {'CO':
                                               {'CH4': 1.0,
                                                'H2': -3.0,
                                                'H2O': 1.0},
                                               'CO2':
                                               {'CH4': 1.0,
                                                'H2': -4.0,
                                                'H2O': 2.0},
                                               'NO':
                                               {'H2': -2.5,
                                                'H2O': 1.0,
                                                'NH3': 1.0},
                                               'O2': {'H2': -2.0, 'H2O': 2.0}},
                     'symbols': ['H', 'C', 'N', 'O']},
                    {'adsorbates': ['NO', 'O2', 'CO'],
                     'references': [('H', 'H2'),
                                    ('C', 'CH4'), ('N', 'NH3'), ('O', 'H2O')],
                     'stoichiometry': [[0.5, 0.0, 0.0, 0.0],
                                       [-2.0, 1.0, 0.0, 0.0],
                                       [-1.5, 0.0, 1.0, 0.0],
                                       [-1.0, 0.0, 0.0, 1.0]],
                     'stoichiometry_factors': {
                        'CO': {'CH4': 1.0, 'H2': -3.0, 'H2O': 1.0},
                        'NO': {'H2': -2.5, 'H2O': 1.0, 'NH3': 1.0},
                        'O2': {'H2': -2.0, 'H2O': 2.0}},
                     'symbols': ['H', 'C', 'N', 'O']},
                    {'adsorbates': ['NO', 'O2', 'CO', 'Cl'],
                     'references': [('H', 'H2'),
                                    ('C', 'CH4'),
                                    ('N', 'NH3'),
                                    ('O', 'H2O'),
                                    ('Cl', 'HCl')],
                     'stoichiometry': [[0.5, 0.0, 0.0, 0.0, 0.0],
                                       [-2.0, 1.0, 0.0, 0.0, 0.0],
                                       [-1.5, 0.0, 1.0, 0.0, 0.0],
                                       [-1.0, 0.0, 0.0, 1.0, 0.0],
                                       [-0.5, 0.0, 0.0, 0.0, 1.0]],
                     'stoichiometry_factors': {
                        'CO': {'CH4': 1.0, 'H2': -3.0, 'H2O': 1.0},
                        'Cl': {'H2': -0.5, 'HCl': 1.0},
                        'NO': {'H2': -2.5, 'H2O': 1.0, 'NH3': 1.0},
                        'O2': {'H2': -2.0, 'H2O': 2.0}},
                     'symbols': ['H', 'C', 'N', 'O', 'Cl']},
                    {'adsorbates': ['NO', 'O2', 'CO', 'Cl', 'HCl', 'Cl2'],
                     'references': [('H', 'H2'),
                                    ('C', 'CH4'),
                                    ('N', 'NH3'),
                                    ('O', 'H2O'),
                                    ('Cl', 'HCl')],
                     'stoichiometry': [[0.5, 0.0, 0.0, 0.0, 0.0],
                                       [-2.0, 1.0, 0.0, 0.0, 0.0],
                                       [-1.5, 0.0, 1.0, 0.0, 0.0],
                                       [-1.0, 0.0, 0.0, 1.0, 0.0],
                                       [-0.5, 0.0, 0.0, 0.0, 1.0]],
                     'stoichiometry_factors': {
                        'CO': {'CH4': 1.0, 'H2': -3.0, 'H2O': 1.0},
                        'Cl': {'H2': -0.5, 'HCl': 1.0},
                        'Cl2': {'H2': -1.0, 'HCl': 2.0},
                        'HCl': {'HCl': 1.0},
                        'NO': {'H2': -2.5, 'H2O': 1.0, 'NH3': 1.0},
                        'O2': {'H2': -2.0, 'H2O': 2.0}},
                     'symbols': ['H', 'C', 'N', 'O', 'Cl']},
                    {'adsorbates': ['NO', 'O2', 'CO', 'Cl', 'HCl', 'CH4'],
                     'references': [('H', 'H2'),
                                    ('C', 'CH4'),
                                    ('N', 'NH3'),
                                    ('O', 'H2O'),
                                    ('Cl', 'HCl')],
                     'stoichiometry': [[0.5, 0.0, 0.0, 0.0, 0.0],
                                       [-2.0, 1.0, 0.0, 0.0, 0.0],
                                       [-1.5, 0.0, 1.0, 0.0, 0.0],
                                       [-1.0, 0.0, 0.0, 1.0, 0.0],
                                       [-0.5, 0.0, 0.0, 0.0, 1.0]],
                     'stoichiometry_factors': {
                        'CH4': {'CH4': 1.0},
                        'CO': {'CH4': 1.0, 'H2': -3.0, 'H2O': 1.0},
                        'Cl': {'H2': -0.5, 'HCl': 1.0},
                        'HCl': {'HCl': 1.0},
                        'NO': {'H2': -2.5, 'H2O': 1.0, 'NH3': 1.0},
                        'O2': {'H2': -2.0, 'H2O': 2.0}},
                     'symbols': ['H', 'C', 'N', 'O', 'Cl']},
                    {'adsorbates':
                        ['NO', 'O2', 'CO', 'CO2', 'Cl', 'HCl', 'CH4', 'H2O'],
                     'references': [('H', 'H2'),
                                    ('C', 'CH4'),
                                    ('N', 'NH3'),
                                    ('O', 'H2O'),
                                    ('Cl', 'HCl')],
                     'stoichiometry': [[0.5, 0.0, 0.0, 0.0, 0.0],
                                       [-2.0, 1.0, 0.0, 0.0, 0.0],
                                       [-1.5, 0.0, 1.0, 0.0, 0.0],
                                       [-1.0, 0.0, 0.0, 1.0, 0.0],
                                       [-0.5, 0.0, 0.0, 0.0, 1.0]],
                     'stoichiometry_factors': {
                            'CH4': {'CH4': 1.0},
                            'CO': {'CH4': 1.0, 'H2': -3.0, 'H2O': 1.0},
                            'CO2': {'CH4': 1.0, 'H2': -4.0, 'H2O': 2.0},
                            'Cl': {'H2': -0.5, 'HCl': 1.0},
                            'H2O': {'H2O': 1.0},
                            'HCl': {'HCl': 1.0},
                            'NO': {'H2': -2.5, 'H2O': 1.0, 'NH3': 1.0},
                            'O2': {'H2': -2.0, 'H2O': 2.0}},
                     'symbols': ['H', 'C', 'N', 'O', 'Cl']},
                    {'adsorbates':
                     ['H3', 'NO', 'O2', 'CO', 'CO2',
                      'Cl', 'HCl', 'CH4', 'H2O'],
                     'references': [('H', 'H2'),
                                    ('C', 'CH4'),
                                    ('N', 'NH3'),
                                    ('O', 'H2O'),
                                    ('Cl', 'HCl')],
                     'stoichiometry': [[0.5, 0.0, 0.0, 0.0, 0.0],
                                       [-2.0, 1.0, 0.0, 0.0, 0.0],
                                       [-1.5, 0.0, 1.0, 0.0, 0.0],
                                       [-1.0, 0.0, 0.0, 1.0, 0.0],
                                       [-0.5, 0.0, 0.0, 0.0, 1.0]],
                     'stoichiometry_factors': {
                         'CH4': {'CH4': 1.0},
                         'CO': {'CH4': 1.0, 'H2': -3.0, 'H2O': 1.0},
                         'CO2': {'CH4': 1.0, 'H2': -4.0, 'H2O': 2.0},
                         'Cl': {'H2': -0.5, 'HCl': 1.0},
                         'H2O': {'H2O': 1.0},
                         'H3': {'H2': 1.5},
                         'HCl': {'HCl': 1.0},
                         'NO': {'H2': -2.5, 'H2O': 1.0, 'NH3': 1.0},
                         'O2': {'H2': -2.0, 'H2O': 2.0}},
                     'symbols': ['H', 'C', 'N', 'O', 'Cl']}]

    parser = optparse.OptionParser()
    parser.add_option('-v', '--verbose', dest='verbose',
                      action='store_true', default=False)
    parser.add_option('-e', '--extra-verbose',
                      dest='extra_verbose', action='store_true', default=False)
    parser.add_option('-n', '--no-test', dest='test',
                      action='store_false', default=True)
    options, args = parser.parse_args()
    verbose = options.verbose
    extra_verbose = options.extra_verbose

    examples = [
        ['NH', 'CO', 'O', 'SH', 'OH', 'CH3'],
        ['CO', 'O'],
        ['CO', 'N'],
        ['NO', 'O2', 'CO', 'CO2'],
        ['NO', 'O2', 'CO'],
        ['NO', 'O2', 'CO', 'Cl'],
        ['NO', 'O2', 'CO', 'Cl', 'HCl', 'Cl2'],
        ['NO', 'O2', 'CO', 'Cl', 'HCl', 'CH4'],
        ['NO', 'O2', 'CO', 'CO2', 'Cl', 'HCl', 'CH4', 'H2O'],
        ['H3', 'NO', 'O2', 'CO', 'CO2', 'Cl', 'HCl', 'CH4', 'H2O'],
    ]

    results = []
    for adsorbates in examples:
        result = {}
        results.append(result)

        if verbose:
            print(
                "\n\n\n\n\n--------------------------------------------------")
            print("ADSORBATES {adsorbates}".format(**locals()))
        result['adsorbates'] = adsorbates

        symbols = molecules2symbols(adsorbates)
        if verbose:
            print("SYMBOlS {symbols}".format(**locals()))
        result['symbols'] = symbols

        references = construct_reference_system(symbols)
        if verbose:
            print("REFERENCES {references}".format(**locals()))
        result['references'] = references

        stoichiometry = get_atomic_stoichiometry(references)
        if verbose:
            print("STOICHIOMETRY MATRIX")
            pprint.pprint(stoichiometry)
        result['stoichiometry'] = stoichiometry

        stoichiometry_factors = get_stoichiometry_factors(
            adsorbates, references)
        if verbose:
            print("STOICHIOMETRY FACTORS")
            pprint.pprint(stoichiometry_factors)
        result['stoichiometry_factors'] = stoichiometry_factors

    if extra_verbose:
        print('\n\n\n\n\n\n')
        pprint.pprint(results)

    if options.test:
        assert results == test_results
