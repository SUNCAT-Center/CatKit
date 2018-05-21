
import os
import sys
import unittest
import tempfile
import pprint
import sqlite3
import json
import click
from click.testing import CliRunner

import cathub


class CommandLineTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_imports(self):
        import cathub
        import cathub.cathubsqlite
        import cathub.postgresql
        import cathub.folderreader
        import cathub.tools
        import cathub.ase_tools
        import cathub.folder2db
        import cathub.db2server
        import cathub.make_folders_template
        import cathub.psql_server_connect
        import cathub.organize

    def test_cli(self):
        runner = CliRunner()
        
        from cathub import reactions, publications        
        runner.invoke(reactions)
        runner.invoke(publications)

        from cathub import make_folders
        runner.invoke(make_folders, ['--create-template', 'template'])
        runner.invoke(make_folders, ['template'])

                      
if __name__ == '__main__':
    unittest.main()
        
