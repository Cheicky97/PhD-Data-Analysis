# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 16:47:52 2024

@author: codiarra
"""
from scipy.constants import e as q, epsilon_0 as eps0, c
from scipy.constants import m_e as me, hbar 
import numpy as np
from rich import print
bohr = 5.2912*1e-11 # m

#
E_homo =2.135 
E_lumo =2.599
tm = 0.16670868
# Energies informations
Eg_exp = 1.84           # in eV (gap exp)
E_00 =  0.68 #0.556     # in eV (position of the max intensity in the spectr)
eps = Eg_exp-E_00


DeltaE = (E_lumo - E_homo + eps) * q # J
E = (E_lumo - E_homo + eps)
dip = tm * bohr**2

#shift flexible 
f12 = (2./3.) * (me / hbar**2) * (DeltaE - eps*q) * dip



# A SECOND EXPRESSION OF A21
A21 = 2 * (DeltaE**3) * (q**2) * dip / (3 * np.pi * hbar**4 *
                                        eps0 * c**3)

tau = 1. / A21

print("[bold green]The exciton lifetime[/]")
print("tau = %.2f [bold red]ps[/]" %(tau * 1e12))
print("\n"*2)