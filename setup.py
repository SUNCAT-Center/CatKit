import setuptools

setuptools.setup(
    name="CatKit",
    version="0.4.2",
    url="https://github.com/SUNCAT-Center/CatKit",

    author="Jacob Boes",
    author_email="jrboes@stanford.edu",

    description="General purpose tools for high-throughput catalysis.",
    license='GPL-3.0',

    packages=[
        'catkit',
        'catkit.pawprint',
        'catgen',
        'catgen.api',
        'catflow'
    ],
    package_dir={'catkit': 'catkit'},
    package_data={
        'catkit': ['data/*.db', 'data/*.json']},

    install_requires=[
        'ase>=3.16.0',
        'matplotlib>=2.2.2',
        'nose>=1.3.7',
        'numpy>=1.14.2',
        'scipy>=1.0.1',
        'networkx>=2.1',
        'sqlalchemy>=1.2.1',
        'spglib>=1.10.3',
        'future>=0.16.0',
        'sklearn',
        'python-coveralls',
    ],
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <4',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
)
