from catkit.gen.route import get_response_reactions
from catkit.gen.route import get_heppel_sellers
import numpy as np

nu = np.array([
    # H2Os, COs, CO2s, H2s, Hs, OHs, Os, HCOOs, H2O, CO, CO2, H2
    [   0,   1,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   0],  # s1
    [   1,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   0,   0],  # s2
    [   0,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   1,   0],  # s3
    [   0,   0,   0,   1,  -2,   0,   0,   0,   0,   0,   0,   0],  # s4
    [   0,   0,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   1],  # s5
    [  -1,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0],  # s6
    [   0,  -1 ,  1,   0,   0,   0,  -1,   0,   0,   0,   0,   0],  # s7
    [   0,  -1,   0,   0,   0,  -1,   0,   1,   0,   0,   0,   0],  # s8
    [   0,   0,   0,   0,   1,  -1,   1,   0,   0,   0,   0,   0],  # s9
    [   0,  -1,   1,   0,   1,  -1,   0,   0,   0,   0,   0,   0],  # s10
    [   0,   0,   1,   0,   1,   0,   0,  -1,   0,   0,   0,   0],  # s11
    [   0,   0,   1,   0,   0,   1,  -1,  -1,   0,   0,   0,   0],  # s12
    [  -1,   0,   0,   1,  -1,   1,   0,   0,   0,   0,   0,   0],  # s14
    [   0,   0,   0,   1,  -1,  -1,   1,   0,   0,   0,   0,   0],  # s15
    [   0,   0,   1,   1,  -1,   0,   0,  -1,   0,   0,   0,   0],  # s17
])

epsilon = np.array([
    # Just a place holder
    [ 0, 0, 0],  # H2OS
    [ 0, 0, 0],  # COS
    [ 0, 0, 0],  # CO2S
    [ 0, 0, 0],  # H2S
    [ 0, 0, 0],  # HS
    [ 0, 0, 0],  # OHS
    [ 0, 0, 0],  # OS
    [ 0, 0, 0],  # HCOOS
    # C, H, O
    [ 0, 2, 1],  # H2O
    [ 1, 0, 1],  # CO
    [ 1, 0, 2],  # CO2
    [ 0, 2, 0],  # H2
])

# Indices of the terminal species
terminal = [8, 9, 10, 11]

RER, species = get_response_reactions(epsilon, terminal, species=True)
sigma = get_heppel_sellers(nu, species[0])

print('Linearly independent set of reaction routes:')
print(sigma, '\n')

print('Overall reaction routes:')
print(np.dot(sigma, nu))
