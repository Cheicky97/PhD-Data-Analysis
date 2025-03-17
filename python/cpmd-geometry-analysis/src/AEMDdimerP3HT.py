# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 10:05:30 2023

@author: codiarra
"""

###############################################################################
# Load and plot data
#from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
#                                  AnnotationBbox)
#from matplotlib.cbook import get_sample_data
#from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import sys
newpath='C:\\Users\\codiarra\\Desktop\\These\\pythonProg'
if newpath not in sys.path:
    sys.path.append(newpath)

import plotCPMD2 as cpmd
###############################################################################
#                           working directory
#path1 = 'Cheick_AEMD\\p3ht\\dodecaRings\\dimer\\3bx31y15A75A\\vib'
#path4 = 'Cheick_AEMD/p3ht/dodecaRings/dimer/3bx31y15A75A/locTher/3'
rac = 'C:\\Users\\codiarra\\Desktop\\These\\Thermique_orga/p3ht/8Rings/dimer/'
path1 = rac + '3bx31y15A75A/1/8AEMD/beta/inplane/phase2/traj1/'
path12= rac + '3bx31y15A75A/1/8AEMD/beta/inplane/phase1/traj3/'
path2 = rac + '3bx31y15A75A/1/8AEMD/beta/interplane/phase2/traj3'


#______________________________________________________________________________
#                               Parameters
#______________________________________________________________________________
N_ions =400     # Tot. nb of ions
N_e = 960       # Tot. nb of e-
N_Alk = 304
N_R = 96
tstep = 4.0    # a.u
cellVol = 31. * 15.75 * 7.6 # Angstrom**3
Lx = 31.0   # Angstrom (inplane)
Lz = 7.6    # Angstrom (interplane)

#______________________________________________________________________________
#                                    AEMD
#______________________________________________________________________________

AEMD = cpmd.AEMD(path=path1, NbIons=N_ions, NbAlk=N_Alk, NbR=N_R, L=Lx,
                 cellVol=cellVol, alkr=True, inplane=True, phase2=True
                 )

AEMD.setMatPltLibDefaults()
AEMD.loadData()
AEMD.plot4totSystm(timeLim=False,endTime=None, delay=False, wait=None,
                  setlimDT=True, limDT=(30,300), xlim=None)
#AEMD.profilT()
#AEMD.mergePhase(path2=path12)
#AEMD.plot4totSystm(timeLim=False, endTime=None, delay=False, wait=None,
#                   setlimDT=False, limDT=(28,300), mergedData=True)
# AEMD.plot4Alkyls(timeLim=False,endTime=None, delay=False, wait=None,
#                    setlimDT=True, limDT=(28,300), xlim=None,
#                    posText=(.0, 50, 12.0))
# AEMD.plot4Rings(timeLim=False,endTime=None, delay=False, wait=None,
#                    setlimDT=True, limDT=(20,300), xlim=None,
#                    posText=(.0, 50, 12.0))