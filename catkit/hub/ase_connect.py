#!/usr/bin/env python
import os
import ase.db


def main():
    host = 'catalysishub.c8gwuc8jwb7l.us-west-2.rds.amazonaws.com'
    server = 'postgresql://catvisitor:' + \
        os.environ['DB_PASSWORD'] + '@{0}:5432/catalysishub'.format(host)
    db = ase.db.connect(server)

    return db
