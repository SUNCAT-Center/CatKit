#!/usr/bin/env python
import re
import os
import json
import requests
import pprint
import ase.db

from catkit.hub.cathubsqlite import CathubSQLite

all_columns = {'reactions': ['chemicalComposition', 'surfaceComposition',
                             'facet', 'sites', 'coverages', 'reactants',
                             'products', 'Equation',
                             'reactionEnergy', 'activationEnergy',
                             'dftCode', 'dftFunctional',
                             'username', 'pubId'],
               'publication': ['pubId', 'title', 'authors', 'journal',
                               'number', 'volume',
                               'pages', 'year', 'publisher', 'doi', 'tags'],
               'publications': ['pubId', 'title', 'authors', 'journal',
                                'volume', 'number',
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

    return execute_graphQL(query_string)


def execute_graphQL(query_string):
    root = 'http://catappdatabase2.herokuapp.com/graphql'
    print('Connecting to database at http://catappdatabase2.herokuapp.com/graphql')
    print('')
    print('Executing query:')
    print('')
    print(query_string)
    print('')
    print('Getting data from server...')
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


def get_reactions(n_results=20, write_db=False, **kwargs):
    """
    Get reactions from server

    Give key value strings as arguments
    """
    queries = {}
    for key, value in kwargs.items():
        key = map_column_names(key)
        if key == 'distinct':
            if value in ['True', 'true']:
                # WARNING: undefined variable name 'query_dict'
                query_dict.update({key: True})
                continue
        try:
            value = int(value)
            queries.update({key: value})
        except BaseException:
            queries.update({key: '{0}'.format(value)})

    subtables = []
    if write_db:
        subtables = ['reactionSystems', 'publication']
    else:
        subtables = []
    data = query(table='reactions', subtables=subtables,
                 columns=all_columns['reactions'],
                 n_results=n_results, queries=queries)

    if not write_db:
        return data

    print('Writing result to Reactions.db')
    for row in data['data']['reactions']['edges']:
        with CathubSQLite('Reactions.db') as db:
            row = row['node']
            key_values = {}
            for key in all_columns['reactions']:
                v = row[key]
                # if isinstance(v, unicode):
                #    v = v.encode('utf-8')
                try:
                    v = json.loads(v)
                except BaseException:
                    pass
                key_values[convert(key)] = v
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

            # publications
            pub_key_values = {}
            row_p = row['publication']
            for key in all_columns['publications']:
                pub_key_values[convert(key)] = row_p[key]
            db.write_publication(pub_key_values)

            # reactions and reaction_systems
            id = db.check(key_values['chemical_composition'],
                          key_values['reaction_energy'])
            if id is None:
                id = db.write(key_values)
            else:
                db.update(id, key_values)
        if ase_ids is not None:
            # Ase structures
            with ase.db.connect('Reactions.db') as ase_db:
                for unique_id in ase_ids.values():
                    if ase_db.count('unique_id={}'.format(unique_id)) == 0:
                        atomsrow = get_atomsrow_by_id(unique_id)
                        ase_db.write(atomsrow)

    return data


def get_publications(**kwargs):
    queries = {}
    for key, value in kwargs.items():
        key = map_column_names(key)
        if key == 'distinct':
            if value in ['True', 'true']:
                # WARNING: undefined variable name 'query_dict'
                query_dict.update({key: True})
                continue
        try:
            value = int(value)
            queries.update({key: value})
        except BaseException:
            queries.update({key: '{0}'.format(value)})

    # WARNING: undefined variable name 'publication_columns'
    return query(table='publications', columns=publication_columns,
                 queries=queries)


def get_ase_db():
    ps = os.environ.get('DB_PASSWORD')
    return ase.db.connect(
        'postgresql://catvisitor:{}@catalysishub.c8gwuc8jwb7l.us-west-2.rds.amazonaws.com:5432/catalysishub'.format(ps))


def get_atomsrow_by_id(unique_id):
    db = get_ase_db()
    row = db.get('unique_id={}'.format(unique_id))
    return row

# def get_atomsrow_by_ids(unique_ids):
#    db = get_ase_db()
#    con = db._connect()
#    cur = con.cursor()


def get_atoms_by_id(unique_id):
    row = get_atomsrow_by_id(unique_id)
    return row.toatoms()


def map_column_names(column):
    mapping = {'surface': 'chemicalComposition'}

    if column in mapping:
        return mapping[column]
    else:
        return column


def convert(name):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


if __name__ == '__main__':
    query = query(table='reactions',
                  columns=['chemicalComposition',
                           'reactants',
                           'products'],
                  n_results=10,
                  queries={'chemicalComposition': "~Pt"})
