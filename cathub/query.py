#!/usr/bin/env python

import requests
import pprint

reaction_columns = ['chemicalComposition', 'surfaceComposition',
                    'facet', 'sites', 'coverages', 'reactants', 'products', 'Equation',
                    'reactionEnergy', 'activationEnergy', 'dftCode', 'dftFunctional',
                    'username', 'pubId']#, 'reactionSystems', 'systems', 'publication']

publication_columns = ['pubId', 'title', 'authors', 'journal', 'year', 'doi', 'tags']


def query(table='reactions',
          columns=['chemicalComposition',
                   'reactants',
                   'products'],
          n_results=10,
          queries={},
          print_output=False):
    
    query_string = graphql_query(table=table,
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

    statement += '    }\n'
    statement += '  }\n'
    statement += '}}'

    return statement


def get_reactions(n_results=20, **kwargs):
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

    return query(table='reactions', columns=reaction_columns,
                 n_results=n_results, queries=queries)

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



