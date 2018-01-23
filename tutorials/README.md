# Generate reactions and molecules

The first step is to create a reaction network. Run `1_generate_gasphase.py` first.

  ```from catkit.pathways import ReactionNetwork
    db_name = 'example-C2H6.db'
    with ReactionNetwork(db_name=db_name) as rn:
        molecules = rn.molecule_search(
            element_pool={'C': 2, 'H': 6},
            multiple_bond_search=False)
        rn.save_molecules(molecules)
  ```

This generates the molecules up to C2H6 and the reaction network connecting these molecules.
Then it stores the molecules and connections in the sqlite3 file `example-C2H6.db`.
You can export the reaction network in catmap notation by running `2_to_catmap.py`.
The code loads `example-C2H6.db` and prints the list of reactions.

  ```from catkit.pathways import catmap_api
     catmap = catmap_api(db_name='example-C2H6.db')
     
     rxn = catmap.rxn_expressions()
  ```

`rxn` is a list of strings, which can be appended to a catmap setup file.
