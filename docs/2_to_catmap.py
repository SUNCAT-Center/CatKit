from catkit.catmap_api import catmap_api

catmap = catmap_api(db_name='example-C2H6.db')

# Get reaction expressions.
rxn = catmap.rxn_expressions()
for r in rxn:
    print(r)
# Get species definitions.
species = catmap.species_definitions()
for s in species.keys():
    print(s, species[s])
