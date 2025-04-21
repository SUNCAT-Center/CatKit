from catkit.gen.route import get_response_reactions
import numpy as np

epsilon = np.array([
    # To keep indexing consistent
    [ 0, 0, 0, 0],  # I1
    [ 0, 0, 0, 0],  # I2
    [ 0, 0, 0, 0],  # I3
    [ 0, 0, 0, 0],  # I4
    [ 0, 0, 0, 0],  # I5
    # C  N  H  O
    [ 1, 0, 4, 0],  # CH4
    [ 0, 1, 0, 1],  # NO
    [ 0, 0, 0, 2],  # O2
    [ 0, 2, 0, 0],  # N2
    [ 1, 0, 0, 1],  # CO
    [ 1, 0, 0, 2],  # CO2
    [ 0, 0, 2, 1],  # H2O
])

terminal = [5, 6, 7, 8, 9, 10, 11]
OR, species = get_response_reactions(epsilon, terminal, species=True)

print('Overall reaction routes:')
print(OR, '\n')

print('Terminal species:')
print(species)
