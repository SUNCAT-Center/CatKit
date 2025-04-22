from catkit.gen.route import get_response_reactions
from catkit.gen.route import get_reaction_routes
from catkit.gen.route import get_heppel_sellers
import numpy as np
np.set_printoptions(threshold=np.inf)

nu = np.array([
    [  1,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0],  # s1
    [  0,  1,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0],  # s2
    [  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  1,  0],  # s3
    [  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0,  0],  # s4
    [  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  1],  # s5
    [ -1,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0],  # s6
    [  0, -1,  1,  0,  0,  0, -1,  0,  0,  0,  0,  0],  # s7
    [  0, -1,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0],  # s8
    [  0,  0,  0,  0,  1, -1,  1,  0,  0,  0,  0,  0],  # s9
    [  0, -1,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0],  # s10
    [  0,  0,  1,  0,  1,  0,  0, -1,  0,  0,  0,  0],  # s11
    [  0,  0,  1,  0,  0,  1, -1, -1,  0,  0,  0,  0],  # s12
    [ -1,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0],  # s13
])

epsilon = np.array([
    # C, H, O
    [ 0, 2, 1],  # SH2O
    [ 1, 0, 1],  # SCO
    [ 1, 0, 2],  # SCO2
    [ 0, 2, 0],  # SH2
    [ 0, 1, 0],  # SH
    [ 0, 1, 1],  # SOH
    [ 0, 0, 1],  # SO
    [ 1, 1, 2],  # SOOCH
    [ 0, 2, 1],  # H2O
    [ 1, 0, 1],  # CO
    [ 1, 0, 2],  # CO2
    [ 0, 2, 0],  # H2
])


# Indices of the species considered terminal
terminal = [8, 9, 10, 11]

RER, species = get_response_reactions(epsilon, terminal, species=True)
sigma = get_heppel_sellers(nu, species[0])
FR, ER = get_reaction_routes(nu, sigma)

print('{} Full reaction routes:'.format(len(FR)))
print(FR, '\n')

print('{} Empty reaction routes:'.format(len(ER)))
print(ER)
