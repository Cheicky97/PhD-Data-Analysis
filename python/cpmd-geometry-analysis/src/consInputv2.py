# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 13:57:59 2023

@author: codiarra
"""
###############################################################################
# Multiple input construction based on geometries in trajectory
#import glob
import os
from shutil import copytree

def recCoord(line, nb, obfile1, obfile2):
    """
    Copy nb following lines (from line) of file 1 into file2.

    Parameters
    ----------
    line : TYPEstr
        line from which the copy will begin (not included).
    nb : int
        number of lines that will be copied.
    obfile1 : _io.TextIOWrapper
        object file 1.
    obfile2 : _io.TextIOWrapper
        object file 2.

    Returns
    -------
    None.

    """
    for j in range(nb):
        line = obfile1.readline()
        obfile2.write(line[5:])
        
def condRec(line, keyword, obfile1, obfile2):
    """
    Copy file1 into file 2 from line to keyword. 

    Parameters
    ----------
    line : str
        line
    keyword : str
       The word that has to be find in the obfile1.
    obfile1 : _io.TextIOWrapper
        The input object file.
    obfile2 : _io.TextIOWrapper
        The output objectfile.

    Returns
    -------
    None.

    """
    while 1:
        line = obfile1.readline()   # first line in input
        obfile2.write(line)         # write it in the created file
        if keyword in line:
            line = obfile1.readline()
            obfile2.write(line)
            break


def createAndfill(nbSpecies, keywords='LMAX=', filename='o-idtbr.in', nbTraj=1,
                  trajFile='TRAJEC_SAMPLED.xyz', source='./0/', target='0',
                  compute='all'):
    """
    Creates a given number directories and CPMD input files for each step in a 
    trajectory (.xyz) file. With the input (resp. directory) base on an input 
    (resp. directory) of reference.

    Parameters
    ----------
    nbSpecies : tup of int
        Each element of the tuple is the number of atoms of given species.
        The order should respect the order in the trajectory file.
    keywords: str
         Corresponds to the line that contains "LMAX=" in CPMD file.
         The default is "LMAX=".
    filename : str, optional
        Name of the CPMD input file that will be used the new input file.
        The default is 'o-idtbr.in'.
    nbTraj : int, optional
        Total number of configuration in trajectroy file. The default is 1.
    trajFile : str, optional
        The name of the trajectory file. The default is 'TRAJEC_SAMPLED.xyz'.
    source : str, optional
        Path to the directory of reference. The default is './0/'.
    target : str, optional
        Name of the target directory (in source directory). The default is '0'.
    compute : str, optional
        decide if and how existing directory will be computed.
        The default is 'all'.
        possible values:
            "all" --> automatically recompute all existing file.
            "partial" --> to manually select the directory to ignore.
            "none" --> ignore all existing file.
            

    Raises
    ------
    FileNotFoundError
        When one file or trajectory in the argument does not exist.

    Returns
    -------
    None.

    """
    path0 = os.getcwd() #
    msg = '{} not found.'
    dirname = os.path.join(source, target) # path to the target directory
    if not os.path.exists(dirname):
        raise FileNotFoundError(msg.format(dirname))
    else:
        if not os.path.exists(os.path.join(dirname, filename)):
            raise FileNotFoundError(msg.format(filename))
    srcFilePath = os.path.join(dirname, filename)
    
    
    if os.path.exists(trajFile):
        trajec = open(os.path.join(path0, trajFile), 'r')
    else:
        raise FileNotFoundError(msg.format(trajFile))
    
    #i = 0
    
    for num in range(1, nbTraj + 1):
        dest = './' + str(num)
        os.chdir(path0)
        if os.path.exists(dest):
            if compute.strip().upper() == 'ALL':
                os.chdir(os.path.join(dest, target))
            elif compute.strip().upper() == 'PARTIAL':
                
                print(dest + " already exist !")
                overwrite = input("Would like to operate on it (yes or no) ? ")
                if overwrite.strip().upper() == 'NO':
                    continue
                else:
                    os.chdir(os.path.join(dest, target))
            elif compute.strip().upper() == "NONE":
                continue
        else:
            copytree(src=source, dst=dest)
            os.chdir(os.path.join(dest, target))
        print(">> The directory", dest[2:], 'has been created.')
        
        lineTraj = trajec.readline() # first line in traj
        srcFile = open(os.path.join(path0, srcFilePath), 'r')
        dstFile = open(filename, 'w')
        #i += 1
        
        if lineTraj:
            lineTraj = trajec.readline()    # second line in traj
            lineIn = srcFile.readline()     # first line in input
            dstFile.write(lineIn)           # write it in the created file
            
            for atom in nbSpecies:
                condRec(lineIn, keywords, srcFile, dstFile)
                recCoord(lineTraj, atom, trajec, dstFile)
            dstFile.write(srcFile.readline())
            srcFile.close()
            dstFile.close()
        else:
            break
        
    trajec.close()
        
        
if __name__ == '__main__':
    #filename = ''
    #         S  O   N   C    H
    atoms = (32, 8, 24, 176, 128)
    
    #createAndfill(atoms, nbTraj=101, compute='None')
    