#!/usr/bin/env python3

import math

# Default values for exfoliated:
# L = 78 # length of the particle (in nm)
# t = 11 # thickness of the particle (in nm)
# E_f = 232 # Young's modulus of the filler (in GPa)
# E_m = 2 # Young's modulus of the matrix (in GPa) to be consistent with FEM

L = 78
t = 11
E_f = 232
nu_f = 0.3
E_m = 2
nu_m = 0.3

def model(fi, L, t, E_f, nu_f, E_m, nu_m):
    alpha = t / L
    g = alpha / (1 - alpha**2)**1.5 * (1 / math.cos(alpha) - alpha * (1 - alpha**2)**0.5)
    S1111 = 1 / 2 / (1 - nu_m) * (1 - 2 * nu_m + (3 * alpha**2 - 1) / (alpha**2 - 1) - (1 - 2 * nu_m + (3 * alpha**2) / (alpha**2 - 1)) * g)
    S2222 = 3 * alpha**2 / 8 / (1 - nu_m) / (alpha**2 - 1) + 1 / 4 / (1 - nu_m) * (1 - 2 * nu_m - 9 / 4 / (alpha**2 - 1)) * g
    S3333 = S2222
    S2233 = 1 / 4 / (1 - nu_m) * (alpha**2 / 2 / (alpha**2 - 1) - (1 - 2 * nu_m + 3 / 4 / (alpha**2 - 1)) * g)
    S3322 = S2233
    S2211 = -alpha**2 / 2 / (1 - nu_m) / (alpha**2 - 1) + g / 4 / (1 - nu_m) * (3 * alpha**2 / (alpha**2 - 1) - (1 - 2 * nu_m))
    S3311 = S2211
    S1122 = -(1 - 2 * nu_m + 1 / (alpha**2 - 1)) / 2 / (1 - nu_m) + (1 - 2 * nu_m + 3 / 2 / (alpha**2 - 1)) * g / 2 / (1 - nu_m)
    S1133 = S1122
    S2323 = (alpha**2 / 2 / (alpha**2 - 1) + (1 - 2 * nu_m - 3 / 4 / (alpha**2 - 1)) * g) / 4 / (1 - nu_m)
    S3232 = S2323
    S1212 = (1 - 2 * nu_m - (alpha**2 + 1) / (alpha**2 - 1) - g / 2 * (1 - 2 * nu_m - 3 * (alpha**2 + 1) / (alpha**2 - 1))) / 4 / (1 - nu_m)
    S1313 = S1212
    
    l_m = nu_m * E_m / (1 + nu_m) / (1 - 2 * nu_m)
    mu_m = E_m / 2 / (1 + nu_m)
    l_f = nu_f * E_f / (1 + nu_f) / (1 - 2 * nu_f)
    mu_f = E_f / 2 / (1 + nu_f)
    D1 = 1 + 2 * (mu_f - mu_m) / (l_f - l_m)
    D2 = (l_m + 2 * mu_m) / (l_f - l_m)
    D3 = l_m / (l_f - l_m)
    B1 = fi * D1 + D2 + (1 - fi) * (D1 * S1111 + 2 * S2211)
    B2 = fi + D3 + (1 - fi) * (D1 * S1122 + S2222 + S2233)
    B3 = fi + D3 + (1 - fi) * (S1111 + (1 + D1) * S2211)
    B4 = fi * D1 + D2 + (1 - fi) * (S1122 + D1 * S2222 + S2233)
    B5 = fi + D3 + (1 - fi) * (S1122 + S2222 + D1 * S2233)
    A1 = D1 * (B4 + B5) - 2 * B2
    A2 = (1 + D1) * B2 - (B4 + B5)
    A3 = B1 - D1 * B3
    A4 = (1 + D1) * B1 - 2 * B3
    A5 = (1 - D1) / (B4 - B5)
    A = 2 * B2 * B3 - B1 * (B4 + B5)
    
    E11_ratio = 1 / (1 + fi * (A1 + 2 * nu_m * A2) / A)
    E22_ratio = 1 / (1 + fi * (-2 * nu_m * A3 + (1 - nu_m) * A4 + (1 + nu_m) * A5 * A) / 2 / A)
    #print('S1111 = ', S1111)
    #print('S2222 = ', S2222)
    #print('S2233 = ', S2233)
    #print('S2211 = ', S2211)
    #print('S1122 = ', S1122)
    #print('S2323 = ', S2323)
    #print('S1212 = ', S1212)
    #print('E11 = ', E11_ratio * E_m)
    #print('E22 = ', E22_ratio * E_m)
    
    return [E11_ratio, E22_ratio]


def getDependenceFi():
    for fi in range(301):
        [E11_ratio, E22_ratio] = model(fi / 10000, L, t, E_f, nu_f, E_m, nu_m)
        print(fi / 10000, E11_ratio * E_m, E22_ratio * E_m)


getDependenceFi()
