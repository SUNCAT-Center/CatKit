<div class="HTML">
<a href='https://travis-ci.org/SUNCAT-Center/CatKit.svg?branch=master'><img src='https://travis-ci.org/SUNCAT-Center/CatKit.svg?branch=master'/></a>
<a href='https://coveralls.io/github/SUNCAT-Center/CatKit?branch=master'><img src='https://coveralls.io/repos/github/SUNCAT-Center/CatKit/badge.svg?branch=master' alt='Coverage Status' /></a>

</div>

Welcome to CatKit! A staging ground for computational tools which are generally useful for catalysis research. The goal of the project is to provide a communal location for those interested in hosting such tools under a common banner. In doing so, we hope to provide the infrastructure to produce more advanced functionality based on modular components of individual contributors.


# Graph-atoms (Gratoms) object

The `Gratoms` object is a child class of the ASE `Atoms` object class. Effectively, this means that it has all the same functionality as an `Atoms` object with some additional features as well.


## Graph information

The `Gratoms` object primary function is to conveniently store and manage graph information.

    from catkit import Gratoms
    from catgen.molecules import get_topologies
    
    atoms0 = get_topologies('C4', saturate=True)[0]
    
    print(type(atoms0))
    print('\nAtomic numbers:\n{}\n'.format(atoms0.numbers))
    print('Connectivity matrix:\n{}'.format(atoms0.connectivity))
    
    atoms1 = atoms0.copy()
    del atoms1[0]
    
    print('\nAtomic numbers:\n{}\n'.format(atoms1.numbers))
    print('Connectivity matrix:\n{}'.format(atoms1.connectivity))
    print('\nAre isomorphs: {}'.format(atoms0.is_isomorph(atoms1)))
    
    atoms2 = atoms1 + Gratoms('C')
    
    print('\nAtomic numbers:\n{}\n'.format(atoms2.numbers))
    print('Connectivity matrix:\n{}'.format(atoms2.connectivity))
    
    atoms2.graph.add_edges_from([(0, 13), (1, 13), (2, 13), (3, 13)])
    
    print('\nAtomic numbers:\n{}\n'.format(atoms2.numbers))
    print('Connectivity matrix:\n{}'.format(atoms2.connectivity))
    
    print('\nAre isomorphs: {}'.format(atoms0.is_isomorph(atoms2)))

<class 'catkit.gratoms.Gratoms'>

Atomic numbers: <br />
[6 6 6 6 1 1 1 1 1 1 1 1 1 1] <br />

Connectivity matrix: <br />
[[0 1 1 1 1 0 0 0 0 0 0 0 0 0] <br />
 [1 0 0 0 0 1 1 1 0 0 0 0 0 0] <br />
 [1 0 0 0 0 0 0 0 1 1 1 0 0 0] <br />
 [1 0 0 0 0 0 0 0 0 0 0 1 1 1] <br />
 [1 0 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 0 1 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 0 1 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 0 1 0 0 0 0 0 0 0 0 0 0]]

Atomic numbers: <br />
[6 6 6 1 1 1 1 1 1 1 1 1 1]

Connectivity matrix: <br />
[[0 0 0 0 1 1 1 0 0 0 0 0 0] <br />
 [0 0 0 0 0 0 0 1 1 1 0 0 0] <br />
 [0 0 0 0 0 0 0 0 0 0 1 1 1] <br />
 [0 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0]]

Are isomorphs: False

Atomic numbers: <br />
[6 6 6 1 1 1 1 1 1 1 1 1 1 6]

Connectivity matrix: <br />
[[0 0 0 0 1 1 1 0 0 0 0 0 0 0] <br />
 [0 0 0 0 0 0 0 1 1 1 0 0 0 0] <br />
 [0 0 0 0 0 0 0 0 0 0 1 1 1 0] <br />
 [0 0 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [1 0 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [1 0 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [1 0 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 0 0 0 0 0 0 0 0 0 0 0 0]]

Atomic numbers: <br />
[6 6 6 1 1 1 1 1 1 1 1 1 1 6]

Connectivity matrix: <br />
[[0 0 0 0 1 1 1 0 0 0 0 0 0 1] <br />
 [0 0 0 0 0 0 0 1 1 1 0 0 0 1] <br />
 [0 0 0 0 0 0 0 0 0 0 1 1 1 1] <br />
 [0 0 0 0 0 0 0 0 0 0 0 0 0 1] <br />
 [1 0 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [1 0 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [1 0 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 1 0 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [0 0 1 0 0 0 0 0 0 0 0 0 0 0] <br />
 [1 1 1 1 0 0 0 0 0 0 0 0 0 0]]

Are isomorphs: True


# CatGen: Catalysis Generator

CatGen is an enumeration module designed to construct various catalytic structures.

-   [X] Gas phase molecules
-   [ ] Bulk structures
-   [X] Surfaces structures
-   [X] Adsorption sites
-   [X] Catalytic structures

It also has functionality for enumeration of other systems relevant to the field of catalysis.

-   [X] Reaction mechanisms
-   [X] Reaction routes

For additional details regarding how the generator operates, including example usage, see the [CatGen readme](./catgen/readme.md).


# Dependencies

CatKit attempts to make use of basic functionalities implemented by existing softwares when possible to extend its capabilities.

A full list of required packaged can be found in [the requirements](./requirements.txt).
