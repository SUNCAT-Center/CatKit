from .postgresql import CathubPostgreSQL
from sys import argv
import os


def main(user):
    db = CathubPostgreSQL(password=os.environ['DB_PASSWORD0'])
    db.create_user(user)


if __name__ == '__main__':
    user = argv[1]
    main(user)
