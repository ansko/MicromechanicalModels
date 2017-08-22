#!/usr/bin/env python3

# Default values for intercalated:
# L = 78 # length of the particle (in nm)
# t = 11 # thickness of the particle (in nm)
# E_f = 232 # Young's modulus of the filler (in GPa)
# E_m = 2 # Young's modulus of the matrix (in GPa) to be consistent with FEM

L = 78
t = 11
E_f = 232
E_m = 2

def model(fi, L, t, E_f, E_m):
    ksi = 2 * L / t
    eta = (E_f / E_m - 1) / (E_f / E_m + ksi)
    ratio = (1 + ksi * eta * fi) / (1 - eta * fi)
    return ratio


def getDependenceFi():
    for fi in range(301):
        ratio = model(fi / 10000, L, t, E_f, E_m)
        print(fi / 10000, ratio * E_m)


getDependenceFi()
