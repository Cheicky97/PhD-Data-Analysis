# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 09:15:34 2023

@author: codiarra
"""

import numpy as np
import matplotlib.pyplot as plt
import os


def gauss(x, mu, height, sigma):
    return height * (1. / (sigma * np.sqrt(2*np.pi))) * np.exp(
                                                -1 * ((x - mu) / sigma)**2 )

def parserTM():
    path = os.getcwd()
    print(path)
    dirn = []
    msg = '{}\t{}\t{}\n'
    
    # loop on the name of directories and files present in the path= './''
    for dirname, _, filenames in os.walk('./'):
        
        # ignore dir that names end in '\\0' or '\\1' (for windows)
        # or '\0' and '\1' (for linux)
        if os.path.split(dirname) in ['0', '1'] in dirname:
            continue
        else:
            # remaining dir names (in the present case) end by '\\2'
            if(dirname[-1] in '2' and dirname[-2] in '\\'):
                dirn.append(dirname) # selection of these dir
                
    with open('transMoment.dat','w') as file:
        for dirname in dirn:
            #os.path.join(path, dirname) contatenate path and dirname
            os.chdir(os.path.join(path, dirname))
            try:
                # first, test if MATRIX.DAT exist in current dir
                f = open("MATRIX.DAT", 'r')
            except OSError:
                # if not, close, print message & continue on followin dir
                f.close()
                print('No MATRIX.DAT in {}'.format(dirname))
                continue
            else:
                # if no exception (i.e MATRIX.DAT exist in current dir)
                with open("MATRIX.DAT", 'r') as matrix:
                    while 1:
                        line = matrix.readline()
                        if line: # True if line != '' (line not empty)
                            line2 = line
                        else:
                            matrix.close()
                            line2 = line2.split()
                            # write energies and trans mom in new file
                            file.write(msg.format(line2[5], line2[6],
                                                  line2[10]))
                            
                            break
if __name__ == '__main__':
    try:
        open('transMoment.dat','r')
    except OSError:
        parserTM()
    else:
        print("File 'transMoment.dat' already exists.")
        decision = 'no' #input("Would you like to update it  (yes or no) ? ")
        if decision.strip().upper() == 'YES':
            parserTM()
            print("The file has been updated.")
        else:
            print("This file will be used to compute lifetime.")

    bohr = 5.2912*1e-11 # m
    
    q = 1.6022e-19 # C
    hbar = 1.055e-34 # J.s
    eps0 = 8.854e-12 #  (c/s)**2 s4 kg−1 m−3
    c = 299792458 # m/s
    me = 9.11e-31 # kg
    D = 8.0e-3 * 1e14 #cm2/s
    
    conv_kin =  1.380649e-23 / 1.6e-19 # kB / eV (eV/K) conv from K to eV
    T = 300 # K
    
    #minE, maxE = 1.6, 2.0
    minE, maxE = 0.2, .7
    npoints = 1000
    
    E_homo, E_lumo, tm = np.loadtxt('transMoment.dat', unpack=True)
    
    Eg_exp = 1.9 # eV (gap exp)
    #Eg_cal = 0.95 # eV (gap DFT)
    E_00 = 0.556 # eV (position of the max intensity in the spectr)
    err = 0 #Eg_exp - E_00
    
    DeltaE = (E_lumo - E_homo + err) * 1.6e-19 # J
    E = (E_lumo - E_homo + err)
    nu = DeltaE / (2 * np.pi * hbar)
    dip = tm * bohr**2
    
    f12 = (2./3.) * (me / hbar**2) * (DeltaE - err* 1.6e-19) * dip

    #print("f12 %.3f" %f)

    
    # A SECOND EXPRESSION OF A21
    A21 = 2 * (DeltaE**3) * (q**2) * dip / (3 * np.pi * hbar**4 * eps0 * c**3)
    
    meanA = np.mean(A21)
    tau = 1. / meanA
    Amax = max(A21) 
    Amin = min(A21)
    tauMax = 1. / Amin
    tauMin = 1. / Amax
    invTau = sum(A21) / len(A21)
    _tau = 1. / A21
    meanTau = np.mean(_tau)
    print(meanA)
    print("tau = %.2f ps" %(tau * 1e12))
    print('L = {:.2f} nm'.format(np.sqrt(D * tau)))
    
    
    # gaussian broadening of the spectra
    g = []
    a = 0
    w = 0.025 # conv_kin * T
    Enrg = np.linspace(minE, maxE, npoints)
    for i in Enrg:
        a = 0
        for j in range(len(E)):
            a += gauss(i, E[j]-err, f12[j], sigma=w)
        g.append(a)
    
    gmax = max(g)


    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(4., 3.), dpi=200)
    ax.plot(Enrg+err, g, alpha=0.8, linewidth=2.0,
               label='gaussian width = %.3f eV' %w)
    ax.set_xlabel('$\Delta E$ (eV)')
    ax.set_ylabel('a. u.')
    #plt.title("Without correction")
    plt.legend(loc=3)
    plt.show()
    plt.savefig('spectr.png')

    
