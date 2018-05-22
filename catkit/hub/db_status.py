import ase
from ase.db import *
from sys import argv
from catkit.hub.catappsqlite import *
import os

catbase = os.environ['data'] + 'winther/'

db = ase.db.connect(catbase + 'atoms.db')
n = db.count('id>0')
print('ASE atoms: ', n)


catapp = CatappSQLite(catbase + 'catapp.db')
con = catapp._connect()
cur = con.cursor()
n = catapp.get_last_id(cur)


print('Catapp:', n)
