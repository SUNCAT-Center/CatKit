from setuptools import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()

with open('readme.org') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='CatKit',
    version='0.3.0',
    packages=['catkit', 'catgen', 'catgen.api'],
    author='Jacob Boes',
    author_email='jrboes@stanford.edu',
    license=license,
    description='Catalysis Automation Tool Kit',
    long_description=readme,
    data_files=['requirements.txt', 'LICENSE'],
    setup_requires=['nose>=1.0'],
    install_requires=required,
)
