# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 11:35:56 2023

@author: codiarra
"""
from math import pi

vol = 28720.91162 #bohr
angstrom = 1e-8 #cm
bohr = 0.529177*angstrom #angstrom
ev = 1.60218e-19 #J
ke = 8.9875e9
hbar = 1.054e-34 #J.s
e = 1.6e-19 # C
E = 1*ev #J
m_e = 9.1e-31 #kg
E_f = [-1.629, -1.029, -0.405]
E_i = -3.334
mu_sq = [0.13572173E+04, 0.25463095E+02, 0.25860629E+02]
omega = [13751.7, 18591.1, 23624] #♣cm-1
c = 299792458 # m/s
moment = 0.24300985E+03*bohr**2 #float(input('moment : '))*bohr
#A21 = (4./3.)*ke*(omega**3/(hbar*c**3))*moment*bohr**2

f12 = (2./3.)*(m_e/(hbar**2))*((E_i-E_f[0])*ev)*mu_sq[0]*(bohr**2)
f13 = (2./3.)*(m_e/(hbar**2))*((E_i-E_f[1])*ev)*mu_sq[1]*(bohr**2)
f14 = (2./3.)*(m_e/(hbar**2))*((E_i-E_f[2])*ev)*mu_sq[2]*(bohr**2)

print ('f1 =', f12/(f12+f13+f14))
print('f2 =', f13/(f12+f13+f14))
print('f3 =', f14/(f12+f13+f14))

# deuxième formulation

f = ((8*pi*m_e*omega[0])/(3*hbar*2*pi))*(mu_sq[0]*(bohr**2))
print(f)


a = mu_sq[0]*(vol**(2/3))

mu1 = 13751.7
mu_sqcm1 =  mu_sq[0]*(bohr**2)

f = mu_sqcm1*(8*pi*m_e*mu1)/(3*hbar*2*pi)