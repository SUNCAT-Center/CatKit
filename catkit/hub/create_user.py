from sys import argv
from catkit.hub.postgresql import CathubPostgreSQL
import os


def main(user):
    db = CathubPostgreSQL(password=os.environ['DB_PASSWORD0'])
    db.create_user(user)


if __name__ == '__main__':
    user = argv[1]
    main(user)
