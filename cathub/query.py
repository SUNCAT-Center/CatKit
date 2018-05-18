#!/usr/bin/env python

import requests
import pprint
from cathubsqlite import CathubSQLite

all_columns = {'reactions': ['chemicalComposition', 'surfaceComposition',
                             'facet', 'sites', 'coverages', 'reactants', 'products', 'Equation',
                             'reactionEnergy', 'activationEnergy', 'dftCode', 'dftFunctional',
                             'username', 'pubId'],               
               'publication': ['pubId', 'title', 'authors', 'journal', 'number', 'volume',
                               'pages',
                               'year', 'publisher', 'doi', 'tags'],
               'publications': ['pubId', 'title', 'authors', 'journal', 'volume', 'number',
                                'pages', 'year', 'publisher', 'doi', 'tags'],
               'reactionSystems': ['name', 'energyCorrection', 'aseId'],
               'publicationSystems': ['pubId', 'aseId']}

def query(table='reactions',
          columns=['chemicalComposition',
                   'reactants',
                   'products'],
          subtables=[],
          n_results=10,
          queries={},
          print_output=False):
    
    query_string = graphql_query(table=table,
                                 subtables=subtables,
                                 columns=columns,
                                 n_results=n_results,
                                 queries=queries)
   # print query_string
    
    return execute_graphQL(query_string)


def execute_graphQL(query_string):
    root = 'http://catappdatabase2.herokuapp.com/graphql'
    print('Connecting to database at http://catappdatabase2.herokuapp.com/graphql')
    print('')
    print('Executing query:')
    print('')
    print(query_string)
    print('')
    data = requests.post(root, {'query': query_string}).json()
    print('Result:')
    print('')
    pprint.pprint(data)
    return data

def graphql_query(table='reactions',
                  subtables=[],
                  columns=['chemicalComposition',
                           'reactants',
                           'products'],
                  n_results=10,
                  queries={}):


    statement = '{'
    statement += '{}(first: {}'.format(table, n_results)
    for key, value in queries.iteritems():
        if isinstance(value, str):
            statement += ', {}: "{}"'.format(key, value)
        elif isinstance(value, bool):
            if value:
                statement += ', {}: true'.format(key)
            else:
                statement += ', {}: false'.format(key)
        else:
            statement += ', {}: {}'.format(key, value)
            
    statement += ') {\n'
    statement += ' totalCount\n  edges {\n    node { \n'
    for column in columns:
        column = map_column_names(column)
        statement += '      {}\n'.format(column)
    for subtable in subtables:
        statement += '      {}'.format(subtable)
        statement += '{\n'
        for column in all_columns[subtable]:
            statement += '        {}\n'.format(column)
        statement += '      }\n'
    statement += '    }\n'
    statement += '  }\n'
    statement += '}}'

    return statement


def get_reactions(n_results=20, write_local=False, **kwargs):
    queries = {}
    for key, value in kwargs.iteritems():
        key = map_column_names(key)
        if key == 'distinct':
            if value in ['True', 'true']:
                query_dict.update({key: True})
                continue
        try:
            value = int(value)
            queries.update({key: value})
        except:
            queries.update({key: '{0}'.format(value)})

    subtables = []
    #if write_local:
    subtables = ['reactionSystems', 'publication']
        
    data = query(table='reactions', subtables=subtables,
                 columns=all_columns['reactions'],
                 n_results=n_results, queries=queries)
    import json
    with CathubSQLite('Reaction.db') as db:
        for row in data['data']['reactions']['edges']:
            row = row['node']
            key_values = {}
            #key_values = {'chemical_composition': row['chemicalComposition'],
            #              'surface_composition': row['surfaceComposition'],
            #              'facet': row['facet'],
            #              'sites': row['sites'],
            #              'coverages': row['coverages'],
            #3              'reactants': row['reactants'],
            #              'products': row['products'],
            #              'reaction_energy': row['reactionEnergy'],
            #              'activation_energy': row['activationEnergy'],
            #              'dft_code': row['dftCode'],
            #              'dft_functional': row['dftFunctional'],
            #              'username': row['username'],
            #              'pub_id': row['pubId']
            #}
            print 'got data'
            for key in all_columns['reactions']:
                v = row[key]
                if isinstance(v, unicode):
                    v = v.encode('utf-8')
                try:
                    v = json.loads(v)
                except:
                    pass
                key_values[convert(key)] = v
            print 'encoded!'
            ase_ids = {}
            energy_corrections = {}
            
            for row_rs in row['reactionSystems']:
                if row_rs['name'] == 'N/A':
                    continue
                ase_ids[row_rs['name']] = row_rs['aseId']
                energy_corrections[row_rs['name']] = row_rs['energyCorrection']

            if not ase_ids:
                ase_ids = None
                energy_corrections = None
            key_values['ase_ids'] = ase_ids
            key_values['energy_corrections'] = ase_ids
            
            id = db.check(key_values['chemical_composition'], key_values['reaction_energy'])
            if id is None:
                id = db.write(key_values)
            else:
                db.update(id, key_values)
            print 'wrote reaction!'

            pub_key_values = {}
            row_p = row['publication']
            for key in all_columns['publications']:
                pub_key_values[convert(key)] = row_p[key]
            db.write_publication(pub_key_values)
            print 'wrote publication!'

    
def get_publications(**kwargs):
    queries = {}
    for key, value in kwargs.iteritems():
        key = map_column_names(key)
        if key == 'distinct':
            if value in ['True', 'true']:
                query_dict.update({key: True})
                continue
        try:
            value = int(value)
            queries.update({key: value})
        except:
            queries.update({key: '{0}'.format(value)})

    return query(table='publications', columns=publication_columns,
                  queries=queries)

def get_structure_by_id(unique_id):
    db = ase.db.connect('postgresql://catvisitor:$DB_PASSWORD@catalysishub.c8gwuc8jwb7l.us-west-2.rds.amazonaws.com:5432/catalysishub')
    
    row = db.get('unique_id={}'.format(unique_id))    
    return row

def get_atoms_by_id(unique_id):
    row = get_structure_by_id(unique_id)
    return row.toatoms()

def map_column_names(column):
    mapping = {'surface': 'chemicalComposition'}

    if column in mapping:
        return mapping[column]
    else:
        return column
    

if __name__ == '__main__':
    query = query(table='reactions',
                 columns=['chemicalComposition',
                          'reactants',
                          'products'],
                 n_results=10,
                 queries={'chemicalComposition': "~Pt"})



def convert(name):
    import re
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()
