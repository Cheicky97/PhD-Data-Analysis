# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 08:55:26 2023

@author: codiarra
"""

#from os import chdir
import numpy as np
#path = 'C://'+''.replace('/', '//')
#chdir(path)
def reg(data, file_geo):
    """
    record geometry

    Parameters
    ----------
    data : 
        DESCRIPTION.
    fileGeo : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    of_g = open(file_geo, 'w')
    msg = '   {:10.6f}      {:10.6f}      {:10.6f}\n'
    itr = 0
    for atn in s:
        for j in range(atn):
            j = j + itr
            of_g.write(msg.format(data[0][j], data[1][j], data[2][j]))
        of_g.write('\n')
        itr += atn
    of_g.close()

*geo, vx, vy, vz = np.loadtxt('GEOMETRY.xyz', unpack=True, skiprows=2,
                              usecols=(1,2,3,4,5,6)
                              )

s = [12, 120, 168]
reg(geo,'geo.txt')
reg((vx, vy, vz), 'velo.txt')
of = open('atm.xyz', 'w')
for i in range(sum(s)):
    of.write(str(i+1)+' ')
of.close()
