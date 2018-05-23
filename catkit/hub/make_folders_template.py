#!/usr/bin/python

import os
import sys
import json

import catkit.hub.ase_tools
from catkit.hub.tools import extract_atoms, check_reaction

username = os.environ['USER']


# ---------publication info------------
def main(
        title='Fancy title',  # work title if not yet published
        authors=['Doe, John', 'Einstein, Albert'],  # name required
        journal='',
        volume='',
        number='',
        pages='',
        year='2017',  # year required
        publisher='',
        doi='',
        tags=[],
        DFT_code='Quantum ESPRESSO',  # for example 'Quantum ESPRESSO'
        DFT_functional='BEEF-vdW',  # For example 'BEEF-vdW'
        #  ---------reaction info-----------
        reactions=[
            {'reactants': ['2.0H2Ogas', '-1.5H2gas', 'star'],
             'products': ['OOHstar@ontop']},
            # {'reactants': ['CCH3'], 'products': ['C', 'CH3']},
            # {'reactants': ['CH3star'], 'products': ['CH3gas', 'star']}
        ],
        energy_corrections={},  # For example: {'H2gas': 0.1},
        bulk_compositions=['Pt', 'Ag'],
        crystal_structures=['fcc', 'hcp'],
        facets=['111'],
        custom_base=None):
    """
    Dear all

    Use this script to make the right structure for your folders.
    Folders will be created automatically when you run the script with python.
    Start by copying the script to a folder in your username,
    and assign the right information to the arguments in the function above.

    You can change the parameters and run the script several times if you,
    for example, are using different functionals or are doing different
    reactions on different surfaces.

    Include the phase if necessary:

    'star' for empty site or adsorbed phase. Only necessary to put 'star' if
    gas phase species are also involved.
    'gas' if in gas phase

    Include the adsorption site if relevant:
    In case of adsorbed species, include the site after 'star' as f.ex
    star@top, star@bridge.

    Remember to include the reactions that gives the adsorption energy of
    reaction intermediates, taking gas phase molecules as references
    (preferably H20, H2, CH4, CO, NH3).

    For example, we can write the adsorption of CH2 as:
          CH4(g) - H2(g) + * -> CH2*

    Here you would have to write ['CH4gas', 'H2gas', 'star'] as 'reactants'
    entry.

    The reaction argument is a list of dictionaries. See examples:

    reactions = [
        {'reactants': ['CH4gas', '-H2gas', 'star'],
         'products': ['CH2star@bridge']},
        {'reactants': ['CH4gas', '-0.5H2gas', 'star'],
         'products': ['CH3star@top']}
        ]

    A new dictionary is required for each reaction, and should include two
    lists, 'reactants' and 'products'. Remember to include a minus sign and
    prefactor in the name when relevant. If your reaction is not balanced,
    you will receive an error when running the script.

    # ---------------surface info---------------------

    facets # If complicated structure: use term you would use in publication
    """

    for reaction in reactions:
        check_reaction(reaction['reactants'], reaction['products'])

    # Set up directories

    if custom_base is not None:
        base = custom_base + '/'
    else:
        catbase = os.path.abspath(os.path.curdir)
        base = '%s/%s/' % (catbase, username)

    if not os.path.exists(base):
        os.mkdir(base)

    publication_shortname = '%s_%s_%s' % (authors[0].split(',')[0].lower(),
                                          title.split()[0].lower(), year)

    publication_base = base + publication_shortname + '/'

    if not os.path.exists(publication_base):
        os.mkdir(publication_base)

    # save publication info to publications.txt
    publication_dict = {'title': title,
                        'authors': authors,
                        'journal': journal,
                        'volume': volume,
                        'number': number,
                        'pages': pages,
                        'year': year,
                        'publisher': publisher,
                        'doi': doi,
                        'tags': tags
                        }

    pub_txt = publication_base + 'publication.txt'
    json.dump(publication_dict, open(pub_txt, 'w'))

    if not len(energy_corrections.keys()) == 0:
        energy_txt = publication_base + 'energy_corrections.txt'
        json.dump(energy_corrections, open(energy_txt, 'wb'))

    def create(path):
        if not os.path.exists(path):
            os.mkdir(path)
        return path

    base = create(publication_base + DFT_code + '/')
    bulk_base = create(base + DFT_functional + '/')
    gas_base = create(bulk_base + '/gas/')

    gas_names = []
    ads_names = []
    from catkit.hub.ase_tools import get_state, clear_state, clear_prefactor
    for i in range(len(reactions)):
        rnames = [r.split('@')[0] for r in reactions[i]['reactants'] +
                  reactions[i]['products']]
        states = [get_state(r) for r in rnames]
        gas_names += [clear_state(clear_prefactor(rnames[i]))
                      for i in range(len(states)) if states[i] == 'gas']

    for name in set(gas_names):
        with open(gas_base + 'MISSING: {}_gas.traj'.format(name), 'w'):
            pass

    for bulk in bulk_compositions:
        for crystal_structure in crystal_structures:
            bulk_name = bulk + '_' + crystal_structure
            facet_base = create(bulk_base + bulk_name + '/')
            with open(facet_base + 'MISSING: {}_bulk.traj'.format(bulk_name),
                      'w'):
                pass

            for facet in facets:
                reaction_base = create(facet_base + facet + '/')
                with open(reaction_base + 'MISSING: empty_slab.traj'
                          .format(bulk_name), 'w'):
                    pass

                for i in range(len(reactions)):
                    rname = '_'.join(reactions[i]['reactants'])
                    pname = '_'.join(reactions[i]['products'])
                    reaction_name = '__'.join([rname, pname])
                    base = create(reaction_base + reaction_name + '/')
                    rnames = [r.split('@')[0] for r in
                              reactions[i]['reactants'] +
                              reactions[i]['products']]
                    states = [get_state(r) for r in rnames]

                    ads_names = [clear_prefactor(clear_state(rnames[i]))
                                 for i in range(len(states))
                                 if states[i] == 'star']

                    for ads in ads_names:
                        if ads == '':
                            continue
                        with open(base + 'MISSING: {}.traj'.format(ads), 'w'):
                            pass


if __name__ == "__main__":
    main()
