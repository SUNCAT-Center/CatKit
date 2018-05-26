import json
from sys import argv

user, pub, DFT, XC, reaction, metal, facet, site, final = argv[1:10]

if site != 'None':
    site = int(site)
else:
    site = None
user_dict = {'user': user,
             'pub_level': int(pub),
             'DFT_level': int(DFT),
             'XC_level': int(XC),
             'reaction_level': int(reaction),
             'metal_level': int(metal),
             'facet_level': int(facet),
             'site_level': site,
             'final_level': int(final)
             }

# WARNING: undefined variable name catbase
user_file = '{0}winther/user_specific/{1}.txt'.format(catbase, user)
json.dump(user_dict, open(user_file, 'w'))
