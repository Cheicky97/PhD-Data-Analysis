# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 14:42:49 2023

@author: codiarra
"""
import numpy as np

def fillBloc(dic):
    """
    Fill the two global blocks for AEMD

    Parameters
    ----------
    dic : dict
        E.g.:
            
            NbS=24 # the first atom specy to occur in the file GEOMETRY.xyz
            NbC=240 # The second one
            NbH=336 # the third

            dic={}
            
            for each specy, we will fix the boundaries of blocks
            for e.g., here the first block include the 1st to the 12th atom.
            While the second from 13th atom to the NbS atom.
            Note : NbS+1 just because np.range(NbS+1) stop at NbS.
            Boundaries are fixed as following:
                  
            dic[NbS] = np.array([1, 12, NbS+1])
            dic[NbC] =  np.array([1, 120, NbC+1]) + NbS
            dic[NbH] =  np.array([1, 168, NbH+1]) + NbC + NbS.

    Returns
    -------
    None.

    """
    bloc1 = []; bloc2 = []
    for keys in dic:
        exch = dic[keys][0]
        for i in range(1, dic[keys].shape[0]):
            min_val = exch
            max_val = dic[keys][i] + i%2 
            exch = max_val
            for val in range(min_val, max_val):
                if i%2 == 1:
                   bloc1.append(val)                 
                elif i%2 == 0:
                    bloc2.append(val)
                
    np.savetxt('BLOC1', bloc1, fmt='%i')
    np.savetxt('BLOC2', bloc2, fmt='%i')
    print('ok')
if __name__ == '__main__':
    NbS = 48
    NbC = 480
    NbH = 672

    dic={}
    dic[NbS] = np.array([1, 16, 41, 48, NbS+1])
    dic[NbC] =  np.array([1, 160, 401, 480, NbC+1]) + NbS
    dic[NbH] =  np.array([1, 224, 561, 672, NbH+1]) + NbC + NbS
    
    fillBloc(dic)    

