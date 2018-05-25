#!/usr/bin/python

import os
import sys
import json
#sys.path.insert(0, "

try:  # sherlock 1 or 2
    sherlock = os.environ['SHERLOCK']
    if sherlock == '1':
        catbase = '/home/winther/data_catapp'
    elif sherlock == '2':
        catbase = '/home/users/winther/data_catapp'
except:  # SUNCAT
    catbase = '/nfs/slac/g/suncatfs/data_catapp'

sys.path.append(catbase)
from tools import extract_atoms, check_reaction


username = os.environ['USER']

# ---------publication info------------


def main(
    title='The Challenge of Electrochemical Ammonia Synthesis: A New Perspective on the Role of Nitrogen Scaling Relations',  # work title if not yet published
    authors=['Montoya, Joseph H.', 'Tsai, Charlie', 'Vojvodic, Aleksandra', 'Norskov, Jens K.'],  # name required
    journal='ChemSusChem',
    volume='8',
    number='13',
    pages='2140-2267',
    year='2015',  # year required
    publisher='',
    doi='10.1002/cssc.201500322',
    tags=[],
    DFT_code='Quantum ESPRESSO',  # for example 'Quantum ESPRESSO'
    DFT_functional='BEEF-vdW',  # For example 'BEEF-vdW'

    #  ---------molecules info-----------

    reactions=[
        {'reactants': ['1.0N2gas', '0.5H2gas', 'star'], 'products': ['NNHstar@ontop']},
        {'reactants': ['1.0N2gas', '1.0H2gas', 'star'], 'products': ['NNH2star@ontop']},
        {'reactants': ['0.5N2gas', 'star'], 'products': ['Nstar@ontop']},
        {'reactants': ['0.5N2gas', '0.5H2gas', 'star'], 'products': ['NHstar@ontop']},
        {'reactants': ['0.5N2gas', '1.0H2gas', 'star'], 'products': ['NH2star@ontop']},
        {'reactants': ['0.5H2gas', 'star'], 'products': ['Hstar@ontop']},

        {'reactants': ['1.0N2gas', '0.5H2gas', 'star'], 'products': ['NNHstar@hollow']},
        {'reactants': ['1.0N2gas', '1.0H2gas', 'star'], 'products': ['NNH2star@hollow']},
        {'reactants': ['0.5N2gas', 'star'], 'products': ['Nstar@hollow']},
        {'reactants': ['0.5N2gas', '0.5H2gas', 'star'], 'products': ['NHstar@hollow']},
        {'reactants': ['0.5N2gas', '1.0H2gas', 'star'], 'products': ['NH2star@hollow']},
        {'reactants': ['0.5H2gas', 'star'], 'products': ['Hstar@hollow']},

        {'reactants': ['1.0N2gas', '0.5H2gas', 'star'], 'products': ['NNHstar@fcc']},
        {'reactants': ['1.0N2gas', '1.0H2gas', 'star'], 'products': ['NNH2star@fcc']},
        {'reactants': ['0.5N2gas', 'star'], 'products': ['Nstar@fcc']},
        {'reactants': ['0.5N2gas', '0.5H2gas', 'star'], 'products': ['NHstar@fcc']},
        {'reactants': ['0.5N2gas', '1.0H2gas', 'star'], 'products': ['NH2star@fcc']},
        {'reactants': ['0.5H2gas', 'star'], 'products': ['Hstar@fcc']},

        {'reactants': ['1.0N2gas', '0.5H2gas', 'star'], 'products': ['NNHstar@bridge']},
        {'reactants': ['1.0N2gas', '1.0H2gas', 'star'], 'products': ['NNH2star@bridge']},
        {'reactants': ['0.5N2gas', 'star'], 'products': ['Nstar@bridge']},
        {'reactants': ['0.5N2gas', '0.5H2gas', 'star'], 'products': ['NHstar@bridge']},
        {'reactants': ['0.5N2gas', '1.0H2gas', 'star'], 'products': ['NH2star@bridge']},
        {'reactants': ['0.5H2gas', 'star'], 'products': ['Hstar@bridge']}

    ],
    energy_corrections={},
    bulk_compositions=['Ag', 'Co', 'Fe', 'Pd', 'Re', 'Ru', 'Au', 'Cu', 'Ir', 'Ni', 'Pt', 'Rh'],
    crystal_structures=['fcc'],
    facets=['111'],
    custom_base=None,
    ):

    """
    Dear all

    Use this script to make the right structure for your folders.
    Folders will be created automatically when you run the script with python.
    Start by copying the script to a folder in your username,
    and assign the right information to the variables below.

    You can change the parameters and run the script several times if you,
    for example, are using different functionals or are doing different reactions
    on different surfaces.


    Include the phase if necessary:

    'star' for empty site or adsorbed phase. Only necessary to put 'star' if
    gas phase species are also involved.
    'gas' if in gas phase

    Include the adsorption site if relevant:
    In case of adsorbed species, include the site after 'star' as f.ex
    star@top, star@bridge.

    Remember to include the adsorption energy of reaction intermediates, taking
    gas phase molecules as references (preferably H20, H2, CH4, CO, NH3).
    For example, we can write the desorption of CH2 as:
    CH2* -> CH4(g) - H2(g) + *
    Here you would have to write 'CH4gas-H2gas' as "products_A" entry.

    See examples:

    reactions = [
        {'reactants': ['CH2star@bridge'], 'products': ['CH4gas', '-H2gas', 'star']},
        {'reactants': ['CH3star@top'], 'products': ['CH4gas', '-0.5H2gas', 'star']}
        ]

    Reaction info is now a list of dictionaries. 
    A new dictionary is required for each reaction, and should include two lists,
    'reactants' and 'products'. Remember to include a minus sign in the name when
    relevant.

    # ---------------surface info---------------------

    facets # If complicated structure: use term you would use in publication
    """

    #  ----------- You're done!------------------------

    # Check reactions
    #assert len(reactants) == len(products_A) == len(products_B)

    for reaction in reactions:
        check_reaction(reaction['reactants'], reaction['products'])

    # Set up directories
    base = '%s/%s/' % (catbase, username)
    if custom_base is not None:
        base = custom_base + '/'

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
    json.dump(publication_dict, open(pub_txt, 'wb'))

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
    from ase_tools import get_state, clear_state, clear_prefactor
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
            with open(facet_base + 'MISSING: {}_bulk.traj'.format(bulk_name), 'w'):
                pass

            for facet in facets:
                reaction_base = create(facet_base + facet + '/')
                with open(reaction_base + 'MISSING: empty_slab.traj'.format(bulk_name), 'w'):
                    pass

                for i in range(len(reactions)):
                    rname = '_'.join(reactions[i]['reactants'])
                    pname = '_'.join(reactions[i]['products'])
                    reaction_name = '__'.join([rname, pname])
                    base = create(reaction_base + reaction_name + '/')      
                    rnames = [r.split('@')[0] for r in reactions[i]['reactants'] + 
                              reactions[i]['products']]
                    states = [get_state(r) for r in rnames]
        
                    ads_names = [clear_prefactor(clear_state(rnames[i])) 
                                 for i in range(len(states)) if states[i] == 'star']

                    for ads in ads_names:
                        if ads == '':
                            continue
                        with open(base + 'MISSING: {}.traj'.format(ads), 'w'):
                            pass

def diagnose(folder):
    miss_list = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if 'MISSING: ' in file:
                traj = file.replace('MISSING: ', '') + '.traj'
                if os.path.isfile(root + '/' + traj):
                    print 'found {}'.format(traj)
                    os.remove(root + '/' + file)
                else:
                    miss_list.append(traj)
    
    if len(miss_list) > 0:
        print 'Files missing'
        print miss_list
    else:
        print 'all files there!'
    return 
    
if __name__ == "__main__":
    main()
