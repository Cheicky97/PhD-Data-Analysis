# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 16:49:17 2023

@author: codiarra
"""

###############################################################################
# To execute before importing lib plotCPMD
import sys
newpath='C:\\Users\\codiarra\\Desktop\\These\\pythonProg'
if newpath not in sys.path:
    sys.path.append(newpath)

import plotCPMD2 as cpmd

syst = 'Sb2S3' #

###############################################################################
#                                   Sb2S3
# NPT 300K 
if syst == 'Sb2S3':
    N_ions = 480
    N_e = 1344
    rac = 'C:\\Users\\codiarra\\Desktop\\These\\Sb2S3/Sb2S3/'
    path2 = rac+'NPT/CMass5E7/stability'
    path1 = rac+'NPT/CMass5E7/'
    path3= rac+'300K.2/28'
    path4= rac+'900K/19'
    file=('out.rdf', 'out2.rdf', 'out3.rdf')
    lab=('Amorphous low T', 'Liquid', 'Amorphous high T')
    sb2s3 = cpmd.StructuralProp(path4)
    sb2s3.setMatPltLibDefaults()
    sb2s3.plotMSDAndDiffCoef()
    #sb2s3.plotPartialgOfR(filenames=file, fileLab=lab, putFilLegend=True)
    #sb2s3.plotBondAngleDist(posAngLab=(0.6, 0.7, 0.10))
    #sb2s3_npt.getInfos()
    #sb2s3_npt.setMatPltLibDefaults()
    #sb2s3_npt.plotStress(repTime=True, tstep=5.0)
    
    #sbNRJ = plotCPMD.ENERGIES(N_ions=N_ions, path=path, N_e=N_e, timeStep=5.0)
    #sbNRJ.setMatPltLibDefaults()
    #sbNRJ.loadData()
    #sbNRJ.repairTime()
    #sbNRJ.plotNVE(NPT=True)