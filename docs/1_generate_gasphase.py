from catkit.pathways import ReactionNetwork

db_name = 'example-C2H6.db'

with ReactionNetwork(db_name=db_name) as rn:

    molecules = rn.molecule_search(
        element_pool={
            'C': 2,
            'H': 6
        }, multiple_bond_search=False)
    rn.save_molecules(molecules)
