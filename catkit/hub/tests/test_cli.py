
import os
import sys
import unittest
import tempfile
import pprint
import sqlite3
import json
import click
from click.testing import CliRunner

import catkit.hub


class CommandLineTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_imports(self):
        import catkit.hub
        import catkit.hub.cathubsqlite
        import catkit.hub.postgresql
        import catkit.hub.folderreader
        import catkit.hub.tools
        import catkit.hub.ase_tools
        import catkit.hub.folder2db
        import catkit.hub.db2server
        import catkit.hub.make_folders_template
        import catkit.hub.psql_server_connect
        import catkit.hub.organize

    def test_cli_query(self):
        runner = CliRunner()
        from catkit.hub.cli import reactions, publications
        runner.invoke(reactions)
        runner.invoke(publications)

    def test_cli_make_folders(self):
        from catkit.hub.cli import make_folders
        runner = CliRunner()
        runner.invoke(make_folders, ['--create-template', 'template'])
        runner.invoke(make_folders, ['template'])

    def test1_cli_read_folders(self):
        from catkit.hub.cli import folder2db
        runner = CliRunner()
        runner.invoke(folder2db, ['aayush/'])

    def test2_cli_db2server(self):
        from catkit.hub.postgresql import CathubPostgreSQL
        from catkit.hub.cli import db2server
        db = CathubPostgreSQL(user='postgres')
        con = db._connect()
        db._initialize(con)
        db.truncate_schema()
        runner = CliRunner()
        runner.invoke(db2server, ['--dbuser=postgres',
                                  'aayush/MontoyaChallenge2015.db'])
    def test3_cli_asedb(self):
        from catkit.hub.cli import ase
        runner = CliRunner()
        runner.invoke(ase, ['--dbuser=postgres', '--dbpassword=None'])

    def test4_show_reactions(self):
        from catkit.hub.cli import show_reactions
        runner = CliRunner()
        runner.invoke(show_reactions, ['aayush/MontoyaChallenge2015.db'])

    def test_reactions(self):
        from catkit.hub.cli import reactions
        runner = CliRunner()
        runner.invoke(reactions, ['-q chemicalComposition=~Co', '-q reactants=H'])

    def test_publications(self):
        from catkit.hub.cli import publications
        runner = CliRunner()
        runner.invoke(publications)

if __name__ == '__main__':
    unittest.main()
