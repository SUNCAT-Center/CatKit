import setuptools

with open('requirements.txt') as f:
    required = f.read().splitlines()

with open('readme.org') as f:
    readme = f.read()

setuptools.setup(
    name="CatKit",
    version="0.3.0",
    url="https://github.com/SUNCAT-Center/CatKit",

    author="Jacob Boes",
    author_email="jrboes@stanford.edu",

    description="General purpose tools for high-throughput catalysis.",
    long_description=readme,

    license='GPL-3.0',

    packages=setuptools.find_packages(),

    install_requires=required,
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <4',

    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Developers',
        'Topic :: High Throughput Catalysis',

        'License :: GPL-3.0',

        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
)
