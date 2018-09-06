import os


def get_pub_id(title, authors, year):
    "construct publication id"
    if len(title.split(' ')) > 1 \
       and title.split(' ')[0].lower() in ['the', 'a']:
        _first_word = title.split(' ')[1].split('_')[0]
    else:
        _first_word = title.split(' ')[0].split('_')[0]

    pub_id = authors[0].split(',')[0].split(' ')[0] + \
                  _first_word + \
                  str(year)
    return pub_id


def extract_atoms(molecule):
    """Return a string with all atoms in molecule"""
    if molecule == '':
        return molecule
    try:
        return float(molecule)
    except BaseException:
        pass
    atoms = ''
    if not molecule[0].isalpha():
        i = 0
        while not molecule[i].isalpha():
            i += 1
        prefactor = float(molecule[:i])
        if prefactor < 0:
            prefactor = abs(prefactor)
            sign = '-'
        else:
            sign = ''
        molecule = molecule[i:]
    else:
        prefactor = 1
        sign = ''
    for k in range(len(molecule)):
        if molecule[k].isdigit():
            for j in range(int(molecule[k]) - 1):
                atoms += molecule[k - 1]
        else:
            atoms += molecule[k]
    if prefactor % 1 == 0:
        atoms *= int(prefactor)
    elif prefactor % 1 == 0.5:
        atoms_sort = sorted(atoms)
        N = len(atoms)
        atoms = ''
        for n in range(N):
            for m in range(int(prefactor - 0.5)):
                atoms += atoms_sort[n]
            if n % 2 == 0:
                atoms += atoms_sort[n]

    return sign + ''.join(sorted(atoms))


def add_atoms(atoms_list):
    add = ''
    sub = ''
    for atoms in atoms_list:
        if isinstance(atoms, float):
            continue
        if len(atoms) > 0 and atoms[0] == '-':
            sub += atoms[1:]
        else:
            add += atoms
    return add.replace(sub, '', 1)


def check_reaction(reactants, products):
    """Check the stoichiometry and format of chemical reaction used for
    folder structure.
    list of reactants -> list of products
    """
    reactant_list = [reactant.split('@')[0].strip(
        'star').strip('gas') for reactant in reactants]
    product_list = [product.split('@')[0].strip(
        'star').strip('gas') for product in products]

    reactant_atoms = [extract_atoms(reactant) for reactant in reactant_list]
    product_atoms = [extract_atoms(product) for product in product_list]

    reactants = add_atoms(reactant_atoms)
    products = add_atoms(product_atoms)

    r_stars = 0
    p_stars = 0

    for i, a in enumerate(reactant_atoms):
        if a == '' or 'star' in reactant_list[i]:
            r_stars += 1
        elif isinstance(a, float):
            r_stars += a
    for a in product_atoms:
        if a == '':
            p_stars += 1
        elif isinstance(a, float):
            p_stars += a
    assert ''.join(sorted(reactants)) == ''.join(sorted(products))


def get_catbase():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    if 'SHERLOCK' in os.environ:
        sherlock = os.environ['SHERLOCK']
        if sherlock == '1':
            catbase = '/home/winther/data_catapp/'
        elif sherlock == '2':
            catbase = '/home/users/winther/data_catapp/'
    elif 'SLAC_ENVIRON' in os.environ:
        catbase = '/nfs/slac/g/suncatfs/data_catapp/'
    else:
        catbase = ''
    return catbase


def get_bases(folder_name):
    user = os.environ['USER']
    data_base = folder_name + '/'
    user_base = folder_name
    return data_base, user, user_base

def clear_state(name):
    name = name.replace('*', '').replace('(g)', '')
    name = name.replace('star', '').replace('gas', '')
    return name


def clear_prefactor(molecule):
    if molecule == '':
        return molecule
    if not molecule[0].isalpha():
        i = 0
        while not molecule[i].isalpha():
            i += 1
        molecule = molecule[i:]
    return molecule


def get_atoms(molecule):
    molecule = clear_state(molecule)
    if molecule == '':
        prefactor = 1
        return molecule, prefactor
    try:
        return '', float(molecule)
    except BaseException:
        pass
    if not molecule[0].isalpha():
        i = 0
        while not molecule[i].isalpha():
            i += 1
        prefactor = molecule[:i]
        if prefactor == '-':
            prefactor = -1
        prefactor = float(prefactor)
        molecule = molecule[i:]
    else:
        prefactor = 1

    temp = ''
    for k in range(len(molecule)):
        if molecule[k].isdigit():
            for j in range(int(molecule[k]) - 1):
                temp += molecule[k - 1]
        else:
            temp += molecule[k]

    molecule = ''.join(sorted(temp))

    return molecule, prefactor


def get_state(name):
    if '*' in name or 'star' in name:
        state = 'star'
    elif 'gas' in name:
        state = 'gas'
    else:
        state = 'star'
    return state
