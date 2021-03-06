#!/usr/bin/env python3

import math

# Default values for exfoliated
# t = 0.7 - thickness of the particle (in nm)
# tau = 5 - thickness of the interphase region (in nm)
# V_d - volume fraction of the filler
# E_f = 232 - Young's modulus of the filler (in GPa)
# E_m = 2 - Young's modulus of the matrix (in GPa) to be consistent with FEM
# E_i = 3 - Young's modulus of the interphase (in GPa) to be consistent with FEM

t = 0.7
tau = 5
V_d = 0.102
E_f = 232
E_m = 2
E_i = 3

def model(t, tau, V_d, E_f, E_m, E_i):
    alpha = math.sqrt((2 * tau / t + 1) * V_d)
    fi = math.sqrt(V_d)
    ratio_reversal = 0
    ratio_reversal += (1 - alpha)
    ratio_reversal += (alpha - fi) / (1 - alpha + alpha * E_i / E_m)
    ratio_reversal += fi / (1 - alpha + (alpha - fi) * E_i / E_m + fi * E_f / E_m)
    return 1 / ratio_reversal

def plotTau():
    for i in range(1, 101):
        print(i, model(t, t * i, V_d, E_f, E_m, E_i) * E_m)

def plotVolumeFraction():
    for i in range(0, 101):
        print(i / 1000, model(t, tau, i / 1000, E_f, E_m, E_i) * E_m)

#plotTau()
plotVolumeFraction()