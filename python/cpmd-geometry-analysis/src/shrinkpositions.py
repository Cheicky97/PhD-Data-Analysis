# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 12:53:18 2024

@author: codiarra
"""

import sys
import numpy as np
newpath='C:\\Users\\codiarra\\Desktop\\These\\pythonProg'
if newpath not in sys.path:
    sys.path.append(newpath)
    
from functions import importData
head = "{:}\n{:}\n"
data = "{:2s}   {:8.4f}    {:8.4f}    {:8.4f}\n"
fileOut ='GEOMETRY_1.xyz'

systm ='GEOMETRY Sb2S3'

axe = 'all' 

formerCelldim = 24.0# in Angstrom
newCelldim = 25.0 # in Angstrom

if type(formerCelldim) is float or int:
    alpha = newCelldim / formerCelldim
else:
    if axe == 'all':
        alpha = []
        for a, b, i in zip(newCelldim, formerCelldim, range(3)):
            alpha[i] = a/b
    elif type(axe) is int:
        alpha = formerCelldim[axe] / newCelldim[axe]

symb, geo = importData(filename='GEOMETRY.xyz', ncols=7, skip=(0,4,5,6),
                       ignore=(0,1), grepColumn=0)


if type(alpha) is float:
    if axe == 'all':
        geoPrime = geo * alpha
    elif type(axe) is int :
        geoPrime = geo * 1.0
        geoPrime[axe] = geoPrime * alpha
elif type(alpha) in (list, tuple):
    geoPrime = np.zeros(geo.shape)
    if axe == 'all':
        for i in range(3):
            geoPrime[i] = geo[i] * alpha[i]
    elif type(axe) is int:
        geoPrime[axe] = geo[axe] * alpha

file = open(fileOut, 'w', encoding='utf-8')
file.write(head.format(len(symb), systm))

for i in range(len(symb)):
    file.write(data.format(symb[i], geoPrime[0, i],
                           geoPrime[1, i], geoPrime[2, i]))
file.close()
    
    