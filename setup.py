from setuptools import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='CatKit',
    version='0.0.1.dev0',
    packages=['catkit'],
    license='GNU General Public License v3.0',
    long_description=open('readme.org').read(),
    install_requires=required,
)
