#!/usr/bin/env python

import requests
import pprint

def main(table, columns, n_results, queries):
    query = graphql_query(table=table,
                          columns=columns,
                          n_results=n_results,
                          queries=queries)

    execute(query)



def execute(query):
    root = 'http://catappdatabase2.herokuapp.com/graphql'
    print('Connecting to database at http://catappdatabase2.herokuapp.com/graphql')
    print('')
    print('Executing query:')
    print('')
    print(query)
    print('')
    data = requests.post(root, {'query': query}).json()
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
    for key, value in queries.items():
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
        statement += '      {}\n'.format(column)

    statement += '    }\n'
    statement += '  }\n'
    statement += '}}'

    return  statement


if __name__ == '__main__':
    query = main(table='reactions',
                 columns=['chemicalComposition',
                          'reactants',
                          'products'],
                 n_results=10,
                 queries={'chemicalComposition': "~Pt"})


