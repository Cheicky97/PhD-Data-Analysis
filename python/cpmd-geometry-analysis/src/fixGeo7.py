# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 09:50:59 2024

@author: codiarra
"""
import numpy as np
from rich import print
from math import radians, degrees

def distance(arr):
    """
    Calculate the distance between

    Parameters
    ----------
    arr : TYPE
        DESCRIPTION.

    Returns
    -------
    dist : TYPE
        DESCRIPTION.

    """
    dist = np.empty((arr.shape[1],))
    for i in range(arr.shape[1]):
        dist[i] = np.linalg.norm(arr[:, i])
    return dist


def transVec(cellVec1, cellVec2, cellVec3):
    trans_vec = np.empty((3,28))
    trans_vec[:, 0] = cellVec1.copy()
    trans_vec[:, 1] = -1*cellVec1
    trans_vec[:, 2] = cellVec2.copy(); trans_vec[:, 3] = -1*cellVec2
    trans_vec[:, 4] = cellVec3.copy(); trans_vec[:, 5] = -1*cellVec3
    trans_vec[:, 6] = cellVec1 + cellVec2; trans_vec[:, 7] = -1*trans_vec[:, 6]
    trans_vec[:, 8] = cellVec1 + cellVec3; trans_vec[:, 9] = -1*trans_vec[:, 8]
    trans_vec[:, 10] = cellVec2 + cellVec3
    trans_vec[:, 11] = -1*trans_vec[:, 10]
    trans_vec[:, 12] = cellVec1 + trans_vec[:, 11]
    trans_vec[:, 13] = -1*trans_vec[:, 12]
    trans_vec[:, 14] = cellVec1 - cellVec2
    trans_vec[:, 15] = -1*trans_vec[:, 14]
    trans_vec[:, 16] = cellVec1 - cellVec3
    trans_vec[:, 17] = -1 * trans_vec[:, 16]
    trans_vec[:, 18] = cellVec2 - cellVec3
    trans_vec[:, 19] = -1 * trans_vec[:, 18]
    trans_vec[:, 20] = cellVec1 + cellVec2 - cellVec3
    trans_vec[:, 21] = -1 * trans_vec[:, 20]
    trans_vec[:, 22] = cellVec1 - cellVec2 + cellVec3
    trans_vec[:, 23] = -1 * trans_vec[:, 22]
    trans_vec[:, 24] = cellVec2 + cellVec3 - cellVec1
    trans_vec[:, 25] = -1 * trans_vec[:, 24]
    trans_vec[:, 26] = cellVec1 + cellVec2 + cellVec3
    trans_vec[:, 27] = -1 * trans_vec[:, 26]
    return trans_vec 


def transGeo(trs_vec, geo):
    geoOut = np.empty((28, *geo.shape))
    for i in range(28):
        geoOut[i] = geo + trs_vec[:, i, np.newaxis]
    return geoOut


def cartToSphr(x, y, z):
    rho = np.sqrt(x*x + y*y)
    ar = np.sqrt(rho**2 + z*z)
    theta = np.arcsin(rho/ar) #np.arctan(rho/2.)
    phi = np.arctan(y/x)
    #print("theta %.2f phi %.2f" %(degrees(theta), degrees(phi)))
    return theta, phi

def cossin(a,b):
    return np.cos(a)*np.sin(b)

def sinsin(a,b):
    return np.sin(a)*np.sin(b)

def coscos(a,b):
    return np.cos(a)*np.cos(b)

def norm(arr, arr2=None):
    return np.sqrt(sum(arr*arr))

def sphericalVec(theta, phi):
    e_r = np.array([cossin(phi, theta), sinsin(theta, phi), np.cos(theta)])
    e_theta = np.array([coscos(phi, theta), cossin(theta, phi),
                        -np.sin(theta)
                        ])
    e_phi = np.array([-np.sin(phi), np.cos(phi), 0.0000])
    return e_r, e_theta, e_phi

def angleUV(u,v):
    """
    Angle en radian entre u et v
    """
    return np.arcsin(norm(np.cross(u, v))/(norm(v) * norm(u)))


class FixGeo():
    """
        Class that try to reconstruct folded geometry.
        That is geometry where where are structure converved atoms and some 
        isolated clusters due to the PBC
    """
    def __init__(self, distEps:dict, geofile:str="GEOMETRY.xyz",  
                 inSet:set = {7, 5, 16, 18,54, 41, 43, 52},
                 vec1:list=[26.0140, 0.0000, 0.0000],
                 vec2:list=[-6.7765, 61.4499, 0.0000],
                 vec3:list=[0.0000, 0.0000, 29.8766],
                 traj:bool=False
                 ):
        """
        

        Parameters
        ----------
        distEps : dict
            The different species and corresponding characteristic interatomic
            distance disreagarding the specy of the with which the bonds are
            made.
        geofile : str, optional
            Name of the geometry file. The default is "GEOMETRY.xyz".
        inSet : set, optional
            Set containing the numbers of some atoms that have not been
            affected by the periodic boundary conditions (PBC).
            Note: for a system consisting of multiple stacks atoms number of
            each stack must be given
            The default is {7, 5, 16, 18,54, 41, 43, 52}.
        vec1 : list, optional
            x axis of the periodic cell as given in CPMD output (in bohr).
            The default is [26.0140, 0.0000, 0.0000].
        vec2 : list, optional
            y-axis of the periodic cell .
            The default is [-6.7765, 61.4499, 0.0000].
        vec3 : list, optional
            z-axis of the periodic cell.
            The default is [0.0000, 0.0000, 29.8766].
        traj : bool, optional
            Indicate if the input file is a trajectory of not.
            The default is False.
            
        Returns
        -------
        None.

        """

        print("[bold red]*[/]"*79)
        print(" "*35+"[bold yellow]Welcome ![/]")
        print("[bold red]*[/]"*79)
        
        self.distEps=distEps
        
        # Read geo file here when the input file is not for a trajectory
        # and initialise some characteristic of the system
        if not traj:
            dtype = np.dtype([("specy", "U1"), ("x", float), ("y", float),
                              ("z", float)])
            print("[magenta]Loading geometry file ...[/]")
            self.specy, *self.geo = np.loadtxt(geofile, unpack=True,
                                               skiprows=2,
                                               usecols=(0, 1, 2, 3),
                                               dtype=dtype
                                               )
            self.geo=np.array(self.geo)
            print("Done.")
            
            # initialisation of the lattice vectors
            Bohr = 0.529177 #Angstrom
            self.cellVec = np.empty((3,3))
            self.cellVec[0] = np.array(vec1) * Bohr #(given in CPMD output)
            self.cellVec[1] = np.array(vec2) * Bohr
            self.cellVec[2] = np.array(vec3) * Bohr
            
            # initialisation of the list containing for each atom,
            # its coordinance i.e. the numbers of first neighbours 
            self.coordAtm = []
            for i in range(len(self.specy)):
                self.coordAtm.append([])
           
            # initialisation of the set of atom that not need to be unfold
            # (inset)
            # first test if a file containing this informations already exists
            try:
                open('inSet.txt', 'r')
            except:
                print("No inSet.txt exists.")
                self.inSet = inSet
            else:
                decision = input("Would you like to use existing inSet.txt ?\n>> ")
                if decision.strip().upper() in ("Y", "YES", "OUI"):
                    self.inSet = set(np.loadtxt('inSet.txt').astype(np.int64))
                else:
                    self.inSet = inSet
            
            # initialisation of the list of the indexes of all atoms    
            self.ind = set(np.arange(len(self.specy)))
            # initialisation of the atom that are not present in inset
            self.outSet = self.ind - self.inSet
            
            # important for saturation of unsat geo
            self.indices=None
            self.tobeSat=None 
        
    def saveGeo(self, geoOut,  _set:set, fileOut:str='GEOMETRY_out.xyz'):
        """
        Save the geometry .xyz consisting only of atoms whose number is
        indicated in _set

        Parameters
        ----------
        geoOut : TYPE
            The geo file from which we desire to construct a partial geo file.
        _set : set
            The set of atom which would be in this file.
        fileOut : str, optional
            The output file. The default is 'GEOMETRY_out.xyz'.

        Returns
        -------
        None.

        """
        head = "{:}\n{:}\n"
        data = "{:2s}   {:8.4f}    {:8.4f}    {:8.4f}\n"
        
        systm =fileOut.split('.')[0]
        
        
        file = open(fileOut, 'w', encoding='utf-8')
        file.write(head.format(len(_set), systm))
        
        for i in list(_set):
            file.write(data.format(self.specy[i], geoOut[0, i],
                                    geoOut[1, i], geoOut[2, i]))
        file.close()
        return print(f"[bold green]{fileOut}[/] has been created.")
    
    def included(self, ind:int, inclus:set, restartGeo:bool):
        """
        Add or not an atom to the ensemble of not dispersed atoms

        Parameters
        ----------
        ind : int
            index .
        inclus : set
            set of not dispersed atoms.
        restartGeo : bool
            indicate if the current calculation has to be done from an updated
            geometry (after some rearrangement)
        Returns
        -------
        inclus : set
            Updated set of not dispersed atoms.

        """
        print("[cyan]Defining sets of dispersed and non dispersed atoms ...[/]")
        
        if restartGeo:
            # geoOut is updated geo file
            geo = self.geoOut
        else:
            geo = self.geo
        # calculate distance of all atoms from the atom n=indexed ind
        dist_ij = distance(geo - geo[:, ind, np.newaxis])
        
        coordAtm = np.where(dist_ij<=self.distEps[self.specy[ind]])
        coordAtm = np.setdiff1d(coordAtm, ind)
        
        # update coordinance information of the given atoms
        self.coordAtm[ind] += list(coordAtm)
        # Check is there has been change in the coordinance of atom ind 
        indInclus = (len(set(coordAtm).intersection(inclus)) != 0)
         
        
        if indInclus:
            # if yes add then the the index of the new neighbor in inclus
            inclus= set(self.coordAtm[ind]).union(inclus)
            
        print("Done.")
        return inclus
       
    def incExc(self, restartGeo:bool=False, MaxIter:int=20, save:bool=True):
        """
        Update list of not dispersed atoms, and save if precised the
        corresponding geometry file 

        Parameters
        ----------
        MaxIter : int, optional
            Maximum iteration. The default is 20.
        save : bool, optional
            If True geometry file will be writted. The default is True.

        Returns
        -------
        None.

        """
        print("[magenta]Updating sets of dispersed and non dispersed atoms ...[/]")
        
        ind=np.arange(len(self.specy))
        # goal is to update inset
        if restartGeo:
            geo = self.geoOut
        else:
            geo = self.geo
            
        for i in range(MaxIter):
            print(f"Iter {i+1}")
            
            for idx in ind:
                self.inSet = self.inSet.union(self.included(idx, self.inSet,
                                                            restartGeo
                                                            )
                                              )
        
        if save:
            # Save the list of not dispersed atoms 
            # (that can be directly used in order to not repeat
            # same operations) 
            np.savetxt('inSet.txt', list(self.inSet))
        
        self.saveGeo(geo, self.inSet, fileOut='GEOMETRY_inSet.xyz')

    def coordAndIndInc(self, ind:int, dist_ijOut:np.array):
        """
        Update coordinance of a given atom

        Parameters
        ----------
        ind : int
            The number of the given atom.
        dist_ijOut : np.array
            distances separating this atom form all the atoms.

        Returns
        -------
        coordAtmIndxOut : np.array
            Atoms that are bonded to atom ind.
        indInclus : bool
            indicate if atom ind make bond with atom in inSet.

        """
        
        coordAtmIndxOut = np.where(dist_ijOut<=self.distEps[self.specy[ind]
                                                            ]
                                   )
        coordAtmIndxOut = np.setdiff1d(coordAtmIndxOut, ind)

        indInclus = (len(set(coordAtmIndxOut).intersection(self.inSet)
                         ) != 0
                     )
        return coordAtmIndxOut, indInclus

    def addInclSuppExc(self, ind:int, geoOut:np.array, bigGeo:np.array):
        """
        update :ist of dispersed (outSet)  and not dispersed atoms (inSet)

        Parameters
        ----------
        ind : int
            indice of the current atom.
        geoOut : np.array
            The updated geometry.
        bigGeo : np.array
            The big geometry containing all the surrounding cells (and their
                                                                   atoms)
            to the main cell.
            Note: just a copies of all 28 possible translation of geoFile
            with respect of cell vectors
        Returns
        -------
        None.

        """
        for i in range(28):
            # calculate distance of atom (ind) all atoms in the ith cell
            dist_ijOut = distance(geoOut - bigGeo[i, :, ind, np.newaxis])
            
            # update the coordinance info, if ind (a dispersed atom) make bond
            # with an inSet atom from the ith cell.
            coordAtmIndxOut, indInclus = self.coordAndIndInc(ind, dist_ijOut)
            
            if indInclus:
                # update inSet by adding ind into not dispersed atoms
                self.inSet.add(ind)
                
                # update inSet
                self.inSet= set(self.coordAtm[ind]).union(self.inSet)
                # replace coord of in in geOut by its position in ith cell
                geoOut[:, ind] = bigGeo[i, :, ind].copy()
                
                # do the same for all the atoms linked to ind and that are
                # inside the main cell 
                for item in list(self.coordAtm[ind]):
                    geoOut[:, item] = bigGeo[i,:, item].copy()
                
                #update the set of dispersed atoms (outSet)
                self.outSet = self.outSet - self.inSet
    
    def suppSetFromInc(self, setToSup:set):
        """
        Delete indicate set from the set of not dispersed atoms 

        Parameters
        ----------
        setToSup : set
            the set to be deleted from inSet.

        Returns
        -------
        None.

        """
        print("Deleteing indicated set atoms from set of non dispersed atoms")
        self.inSet = self.inSet - setToSup
        self.out = self.ind - self.inSet
        print("Job done.")
    
    def reconstruct(self, MaxIter:int=100, saveAll:bool=True,
                    restartGeo=False
                    ):
        """
        Reconstruct the molecule (or cystal)

        Parameters
        ----------
        MaxIter : int, optional
            total number of iterations. The default is 100.
        saveAll : bool, optional
            indicate if the final geometry has to be saved.
            The default is True.
        restartGeo : bool, optional
            Indicate if the calculation must restart from an updated geometry.
            The default is False.

        Returns
        -------
        None.

        """
        print("[bold red]Reconstructing the molecule ...[/]")
        # Calculate Outset
        self.outSet = self.ind - self.inSet
        
        # if restart geo restart from the update geometry (geoOut)
        if not restartGeo:
            self.geoOut = self.geo.copy()
        
        # Calculate all the 28 translation vectors
        self.trs_vec = transVec(self.cellVec[0], self.cellVec[1],
                              self.cellVec[2]
                              )
        
        # calculate bigGeo containg all surrounding cell geometries
        self.bigGeo = transGeo(self.trs_vec, self.geoOut)
        
        for i in range(MaxIter):
            print(f"iteration {i+1}")
            # update inSet and outSet
            for ind in self.outSet:
                self.addInclSuppExc(ind, self.geoOut, self.bigGeo)
            print("len Inset", len(self.inSet), "len outSet", len(self.outSet))
        
        print("Saving inSet ...")
        np.savetxt('inSet.txt', list(self.inSet))
        
        # save the geometry of the inSet atoms
        self.saveGeo(self.geoOut, self.inSet, fileOut='GEOMETRY_inSet2.xyz')
        
        if saveAll:
            self.saveGeo(self.geoOut, self.ind, fileOut='GEOMETRY_out.xyz')
        
    def calcRelGeo(self, save:bool=True):
        """
        calculate the translation vectors need for each atoms in order to
        unfolded the periodic boundary conditions

        Parameters
        ----------
        save : bool, optional
            Save translation vectors for ulterior uses and also to evitate to
            recalculate them. The default is True.

        Returns
        -------
        None.

        """
        # definition of the reference geo
        self.rGeo = self.geoOut - self.geo
        
        if save:
            self.saveGeo(self.rGeo, self.ind, fileOut="GEOMETRY_R.xyz")

    def fixTraj(self,fileTraj:str="TRAJEC.xyz"):
        """
        Suppress the periodic boundary conditions from the trajectory

        Parameters
        ----------
        fileTraj : str, optional
            The trajectory file. The default is "TRAJEC.xyz".

        """
        print("[magenta]Fixing the trajectory ...[/]")
        fileOut= fileTraj.split(".")[0]+"_out.xyz"
        file = open(fileOut, 'w', encoding='utf-8')
        
        data = "{:2s}   {:8.4f}    {:8.4f}    {:8.4f}\n"
        print("Verifying if there is existing GEOMETRY_R.xyz")
        try :
            open("GEOMETRY_R.xyz", "r")
        except:
            rGeo = self.rGeo
        else:
            dtype = np.dtype([("x", float), ("y", float), ("z", float)])
            
            rGeo = np.loadtxt("GEOMETRY_R.xyz", unpack=True, skiprows=2,
                              usecols=(1, 2, 3), dtype=dtype
                              )
            rGeo = np.array(rGeo)
            
        with open(fileTraj, "r") as traj:
            for line in traj:
                line += traj.readline()
                file.write(line)
                for i in range(rGeo.shape[1]):
                    line = traj.readline().strip().split()
                    coord = np.asanyarray(line[1:]).astype("float64")                  
    
                    coord = coord + rGeo[:, i]
                    
                    file.write(data.format(line[0], coord[0], coord[1],
                                            coord[2]
                                            )
                                )
        file.close()
        return print(f"A file [bold green]{fileOut}[/] has been created.")
    
    def findInd(self, geoRef:np.array, eps:float=1):
        indices = []
        for ind in range(geoRef.shape[1]):
            dist_ij = distance(self.geo - geoRef[:, ind, np.newaxis])
            coordAtm = np.where(dist_ij<=eps)

            try:
                coordAtm[0][0]
            except:
                continue
            else:
                indices.append(coordAtm[0][0])
        np.savetxt('indices.txt', indices)
        return indices
    
    
    def unSaturatedGeo(self, ind:int, geo:np.array):
        """
        give geoUn the geometry  of the unsaturaded Molecule, given the index
        of atoms in geo that are also in the unsaturade geometry, for example,
        molecule with delete sidechains

        Parameters
        ----------
        ind : int
            DESCRIPTION.
        geo : np.array
            raw geometry (with side chains for example) 
        Returns
        -------
        None.

        """
        geoUn = np.empty((3,len(ind)))
        for i in range(len(ind)):
            geoUn[:,i] = geo[:,ind[i]]
        return geoUn

    def unsaturatedC(self):
        """
        Search for unsaturated carbon atoms 

        Returns
        -------
        None.

        """
        self.tobeSat = []
        for i in range(len(self.specyUn)):
            if self.specyUn[i] == 'C':
                dist_ij = distance(self.geoUn - self.geoUn[:, i, np.newaxis])
                coordAtm = np.where(dist_ij<=distEps['C'])
                vois=list(coordAtm[0])
                if len(vois) == 2:
                    #print(coordAtm[0])
                    self.tobeSat.append([vois, i, round(dist_ij[vois[1]],2),
                                         self.specyUn[vois[1]]]
                                        )
        of = open("tobeSat", "w")
        of.write(str(self.tobeSat))
        of.close()
           
    
    def addHTo(self, coord, d=1.10,alpha=109.47, beta=109.47,
               e_r:np.array=None, e_theta:np.array=None, e_phi:np.array=None):
        """
        Add an Hydrogen to saturad carbon atoms

        Parameters
        ----------
        coord : TYPE
            DESCRIPTION.
        d : TYPE, optional
            DESCRIPTION. The default is 1.10.
        alpha : TYPE, optional
            DESCRIPTION. The default is 109.47.
        beta : TYPE, optional
            DESCRIPTION. The default is 109.47.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        alpha = radians(alpha) 
        beta = radians(beta)
        e_H = (np.cos(alpha)*e_r + cossin(beta, alpha)*e_theta
               + sinsin(alpha, beta)*e_phi)
        return d * e_H
    
    def satur(self, alpha=109.470, beta=120.000, d=1.10, MaxIter:int=4):
        """
        saturate the unsaturated carbon atoms (only to methyl for the moment)

        Parameters
        ----------
        alpha : floar, optional
            The vertical angle. The default is 109.47.
        beta : TYPE, optional
            Azimuth angle in degree. The default is 120.
        d : float, optional
            The C-H distance. The default is 1.10.
        MaxIter : int, optional
            The total number of iteration for the correction in theta in order 
            to have r and e_r collinear
            The default is 3.
        Returns
        -------
        None.

        """
        nbInsat = len(self.tobeSat)
        nbofmissAtm=3
        self.satH =[]
        #np.random.seed(10)
        for i in range(nbInsat):
            _beta = np.random.randint(0.0, 360)
            atmInfo = self.tobeSat[i]
            # initialisation of the atomic positions r0 and r1
            r0 =  self.geoUn[:,atmInfo[1]]
            r1 = self.geoUn[:, atmInfo[0][1]]
            
            # calculation of azimutal angles
            r=(r1-r0)/norm(r1-r0)
            theta, phi = cartToSphr(*r)
            
            #calculation of the sperical vectors basis
            e_r, e_theta, e_phi = sphericalVec(theta, phi)
            
            # a problem with current script is that the calculated e_r is not
            # // to the vector r so we force it by performing successive
            # correction the angle theta
            
            gamma = angleUV(r, e_r) # angle between r and the vecotr e_r
            
            # correction on the angle theta 
            theta = theta + gamma
            # recalculate basis vectors and gamma 
            e_r, e_theta, e_phi = sphericalVec(theta, phi)
            gamma = angleUV(r, e_r)
            
            for it in range(MaxIter):
                
                if gamma > 1e-4:
                    theta += gamma
                    e_r, e_theta, e_phi = sphericalVec(theta, phi)
                    gamma = angleUV(r, e_r)
                else:
                    break
            
            # sometimes for some pair of carbon atoms e_r and r are
            # antiparallel when alpha=0, so correct such case by doing rotation
            # of pi - alpha insted of directly alpha
            
            if r@e_r < 0:
                _alpha = 180.0000 - alpha
            else:
                _alpha = alpha
            
            
            for j in range(nbofmissAtm):
                H_i = self.addHTo(coord=r, d=d, alpha=_alpha,
                                  beta=beta*j + _beta, e_r=e_r,e_theta=e_theta,
                                  e_phi=e_phi
                                  )
                
                self.satH.append(H_i + r0)
        of=open("lenSatH", "w")
        of.write(str(len(self.satH)))
        of.close()

    def delSideChainsFromGeo(self,fileRef:str=None, fileGeo:str=None,
                             eps:float=1e-3, saveUnSatGeo:bool=False,
                             saveSatGeo:bool=True
                             ):
        """
        Delete side chains form geo

        Parameters
        ----------
        fileRef : str, optional
            Exemple of geofile without sidechains.
        fileGeo : str, optional
            geometry from which we want to suppress sidechains.
        eps : float, optional
            Tolerance for interatomic distances. The default is 1e-3.
        saveUnSatGeo : bool, optional
            indicate if to save geometry without sidechains
            The default is False.
        saveSatGeo : bool, optional,
            If True save staurate carbons then save the final geometry.
        Returns
        -------
        None.

        """
        print("[magenta]Removing side Chains from the Trajectory ...[/]")
        
        # read geo file
        if fileGeo != None:
            
            dtype = np.dtype([("specy", "U1"), ("x", float), ("y", float),
                              ("z", float)])
            print("[magenta]Loading geometry file ...[/]")
            self.specy, *self.geo = np.loadtxt(fileGeo, unpack=True,
                                               skiprows=2,
                                               usecols=(0, 1, 2, 3),
                                               dtype=dtype
                                               )
            self.geo=np.array(self.geo)
        else:
            fileGeo="GEOMETRY.xyz"
    

        # calculating indices 
        if self.indices == None:   
            try :
                open("indices.txt", "r")
            except:
                dtype = np.dtype([("x", float), ("y", float), ("z", float)])
                
                rGeo = np.loadtxt(fileRef, unpack=True, skiprows=2,
                                  usecols=(1, 2, 3), dtype=dtype
                                  )
                rGeo = np.array(rGeo) 
                self.indices = self.findInd(rGeo, eps)
            else:
                decis=input("File indices.txt already exist to want to use it ?\n>> ")
                if decis in ["y", "yes"]:
                    self.indices = np.loadtxt("indices.txt").astype('int64')
                else:
                    dtype = np.dtype([("x", float), ("y", float), ("z", float)])
                    
                    rGeo = np.loadtxt(fileRef, unpack=True, skiprows=2,
                                      usecols=(1, 2, 3), dtype=dtype
                                      )
                    rGeo = np.array(rGeo)
                    self.indices = self.findInd(rGeo, eps)
                
        # calculating unsaturated geometry (without sidechains)
        self.geoUn = self.unSaturatedGeo(self.indices, self.geo)
        self.specyUn =[]
        for i in range(len(self.indices)):
             self.specyUn.append(self.specy[self.indices[i]])
         
        # calculating in formation about unsaturated carbons
        self.unsaturatedC()
        
        if saveUnSatGeo:
            self.saveGeo(geoOut=self.geoUn,
                         _set=set(np.arange(len(self.indices)))
                         )
        # saturation with Hydorgens
        if saveSatGeo:
            self.satur()
            
            data = "{:2s}   {:8.4f}    {:8.4f}    {:8.4f}\n"
            fileOut= fileGeo.split(".")[0]+"_out.xyz"
            file = open(fileOut, 'w', encoding='utf-8')
            lenGeoUn = self.geoUn.shape[1]
            lenSatH = len(self.satH) 
            NbAtm = ( lenSatH + lenGeoUn)
            file.write(f"{NbAtm}\nGEOMETRY WITH METHYLS\n" )
            for i in range(lenGeoUn):
                file.write(data.format(self.specyUn[i], self.geoUn[0, i],
                                       self.geoUn[1, i], self.geoUn[2, i])
                           )
            for i in range(lenSatH):
                a = self.satH[i]
                file.write(data.format("H", a[0], a[1], a[2]))

        
    def delSideChainsFromTraj(self,fileTraj:str):
        data = "{:2s}   {:8.4f}    {:8.4f}    {:8.4f}\n"
        
        fileOut= fileTraj.split(".")[0]+"_out.xyz"
        file = open(fileOut, 'w', encoding='utf-8')
        lenGeoUn = self.geoUn.shape[1]
        lenSatH = len(self.satH) 
        NbAtm = ( lenSatH + lenGeoUn)
        with open(fileTraj, "r") as traj:
            for line in traj:
                x = []
                y = []
                z = []
                line = f"{NbAtm}\n" + traj.readline()
                file.write(line)
                
                for i in range(self.geo.shape[1]):
                    line = traj.readline().strip().split()
                    x.append(float(line[1]))
                    y.append(float(line[2]))
                    z.append(float(line[3]))
    
                geo = np.array([x, y, z])#.astype('float64')
                
                self.geoUn = self.unSaturatedGeo(self.indices, geo)
                #self.geo.unsaturatedC()
                self.satur()
                print(self.specyUn[0])
                for i in range(lenGeoUn):
                    file.write(data.format(self.specyUn[i], self.geoUn[0, i],
                                           self.geoUn[1, i], self.geoUn[2, i])
                               )
                for i in range(lenSatH):
                    a = self.satH[i]
                    file.write(data.format("H", a[0], a[1], a[2]))
        file.close()
                


if __name__=="__main__":
        distEps = {}
        distEps['S'] = 1.95
        distEps['O'] = 1.8
        distEps['N'] = 1.8
        distEps['C'] = 1.9
        distEps['H'] = 1.6
        # geo = FixGeo(geofile="GEOMETRY.xyz", distEps=distEps)
        # geo.delSideChainsFromTraj(fileTraj="TRAJ.xyz")
        # geo.unsaturatedC()
        # geo.satur(alpha=109.47, beta=120)
        # geo.incExc(MaxIter=20, save=False)
        # geo.reconstruct(MaxIter=5, saveAll=False)
        # geo.incExc(restartGeo=True, save=False)
        # geo.reconstruct(MaxIter=3, restartGeo=True, saveAll=False)
        geo = FixGeo(geofile="GEOMETRY0.xyz", distEps=distEps)
        geo.delSideChainsFromGeo(fileRef="GEOMETRY_mod.xyz",
                                  fileGeo="GEOMETRY.xyz")
        

        

        
        
        