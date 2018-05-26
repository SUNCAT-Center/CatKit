#!/usr/bin/env python
import os


def main(user):
    os.system('psql --host=catalysishub.c8gwuc8jwb7l.us-west-2.rds.amazonaws.com --port=5432 --username={} --dbname=catalysishub --password'.format(user))


if __name__ == '__main__':
    from sys import argv
    user = argv[1]
    main(user)
