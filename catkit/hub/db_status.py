from .cathubsqlite import CathubSQLite
import ase
import os

# WARNING: Use of unspecified environment variables
# This module throws exceptions when attempting to produce automatic
# documentation
catbase = os.environ['data'] + 'winther/'

db = ase.db.connect(catbase + 'atoms.db')
n = db.count('id>0')
print('ASE atoms: ', n)


catapp = CathubSQLite(catbase + 'catapp.db')
con = catapp._connect()
cur = con.cursor()
n = catapp.get_last_id(cur)

print('Catapp:', n)
