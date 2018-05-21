import setuptools

with open('requirements.txt', 'r') as f:
    requirements = f.readlines()

with open('readme.org', 'r') as f:
    readme = f.read()

setuptools.setup(
    name="CatKit",
    version="0.4.4",
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
        'catkit.gen.api',
        'catkit.flow',
        'cathub',
        'cathub.ase_tools'
    ],
    package_dir={'catkit': 'catkit'},
    entry_points='''
         [console_scripts]
         cathub=cathub:cli
      ''',
    install_requires=requirements,
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
