import setuptools

with open('requirements.txt', 'r') as f:
    requirements = f.readlines()
    requirements = [f for f in requirements if 'git' not in f]
    requirements += ['ase>=3.16.2']

with open('readme.org', 'r') as f:
    readme = f.read()

setuptools.setup(
    name="CatKit",
    version="0.5.4",
    url="https://github.com/SUNCAT-Center/CatKit",

    author="Jacob Boes",
    author_email="jrboes@stanford.edu",

    description="General purpose tools for high-throughput catalysis.",
    long_description=readme,
    license='GPL-3.0',

    packages=[
        'catkit',
        'catkit.pawprint',
        'catkit.gen',
        'catkit.gen.analysis',
        'catkit.gen.utils',
        'catkit.flow',
        'catkit.hub',
        'catkit.hub.ase_tools'
    ],
    package_dir={'catkit': 'catkit'},
    entry_points='''
         [console_scripts]
         cathub=catkit.hub.cli:cli
      ''',
    install_requires=requirements,
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <4',
    dependency_links=[
        'git+https://gitlab.com/ase/ase.git@d441dd6a1c71a2e1a925043e6974a9b3ae961854#egg=ase-3.16.3'
    ],

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
