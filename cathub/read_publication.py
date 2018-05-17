from sys import argv
from postgresql import CathubPostgreSQL
import os
import json

db = CathubPostgreSQL(password=os.environ['DB_PASSWORD0'])

#n_pub = db.status('publication')

columns, rows = db.read(table='publication', id='all')


for row in rows[:0]:
    #print columns, row
    data = {}
    for i, c in enumerate(columns):
        #print c[0], row[0][i]
        data.update({c[0]: row[i]})
    del data['pubtextsearch']
    print data
    pub_txt = '../publications/{}.txt'.format(data['pub_id'])
    json.dump(data, open(pub_txt, 'wb'))    



#pub_id = 'mgfieldslanders2018'
pub_id = 'MichalLixCoO22017'
pub_txt = '../publications/{}.txt'.format(pub_id)

print pub_txt
data = json.load(open(pub_txt, 'r'))
#data = json.dumps(data)

db.update_publication(data)
