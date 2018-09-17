import numpy as np


def get_modified_spin_symbols(numbers, magmoms):
    """Return a representation of atomic symbols which is
    unique to the magnetic moment as well.

    This is effectivly creating a single integer which contains the
    atomic number and the magnetic moment multiplied by 10.

    Parameters
    ----------
    numbers : ndarray (N,)
        Atomic numbers to be joined with the magnetic moments.
    magmoms : ndarray (N,)
        Magnetic moments to be joined to the atomic numbers.

    Returns
    -------
    spin_mod_symbols : ndarray (N,)
        The spin modified symbols representation for each atom.
    """
    spin_mod_symbols = numbers.copy()
    magmoms = magmoms * 10
    magmoms = magmoms.astype(int)

    sign = np.sign(magmoms)
    spin_mod_symbols *= 1000
    spin_mod_symbols += np.abs(magmoms)
    ind = np.where(sign)
    spin_mod_symbols[ind] *= sign[ind]

    return spin_mod_symbols


def get_unmodified_spin_symbols(spin_mod_symbols):
    """Return the origional atomic numbers and magnetic moments from
    the get_modified_spin_symbols function.

    Parameters
    ----------
    spin_mod_symbols : ndarray (N,)
        Joint symbol and spin representation of an atomic structure.

    Returns
    -------
    symbols : ndarray (N,)
        The origional atomic numbers of the atoms object.
    magmoms : ndarray (N,)
        The magnetic moments of the origional atoms object.
    """
    symbols = spin_mod_symbols.copy()
    sign = np.sign(symbols)

    symbols *= sign
    magmoms = symbols % 1000
    symbols -= magmoms
    magmoms = magmoms * (sign / 10)
    symbols //= 1000

    return symbols, magmoms
