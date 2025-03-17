# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:15:53 2023

@author: codiarra
"""
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from scipy.optimize import curve_fit
from os import chdir, getcwd
from myPlots.functions import*
from myPlots.functions import conv_fic, conv_kin, convKToeV, sliceArray
from myPlots.functions import transGeo, distance, finMaxNrj, reorganize_data
from myPlots.functions import cartToSphr, cossin, sinsin, norm, sphericalVec, angleUV
from myPlots.functions import parserTM, plot_AEMD, matplotlib_params
from math import radians#, degrees
from scipy.constants import h, c, eV

from rich import print 
import pandas as pd


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

def _meanT(T):
    T_out = np.zeros(T.shape)
    for i in range(T.shape[0]):
        T_out[i] = T[0:i+1].mean()
    return T_out

###############################################################################
#                                   ENERGIES                                  #
###############################################################################
class ENERGIES(object):
    """
    Class containing all the important plot to visualize energies data
    """

    def __init__(self, path, N_ions, N_e, timeStep=1.):
        """
        

        Parameters
        ----------
        path : str
            path to the directory of interest.
        N_ions : int
            total number of ions in the system.
        N_e : int
            total number of electron in the system.
        timeStep : float, optional
            timeStep.
            The default is 1.0
        Returns
        -------
        None.

        """
        self.path=path
        chdir(self.path)
        self.tstep = timeStep * 0.024 * 1e-3
        self.NbIons, self.NbElec = N_ions, N_e
        print("[bold magenta][cyan]*[/][/]"*79,
              "[bold magenta][cyan]*[/][/]"+' '*35+
              '[bold magenta]ENERGIES[/]'+
              ' '*34+"[bold magenta][cyan]*[/][/]",
              "[bold magenta][cyan]*[/][/]"*79,
              sep='\n')  
        print(
            """
            For informations about all the available methods :
                -> Apply getInfos()
            For informations about the arguments of a given method :
                -> Apply infoMeth(obj.method); with obj the object you created
                    and its method.
                                        
                    Or simply use :
                        >> print(instance.method.__doc__)
            """
            )
            
    def infoMeth(self, method):
        """
        Print the documenation about a given method.

        Parameters
        ----------
        method : method
            The method of interest.
            e.g: if func() is a method of the class object h / h = Test()
            h.infoMeth(h.func) print the documentation of the method func()
        Returns
        -------
        None.

        """
        print(method.__doc__)
        
    def getInfos(self):
        msg ="""
        This Class Contains the following methods:
            changePath() :
                If you want to change the path from the one you
                defined when creating the object structuralProp().
            
            setMatPltLibDefaults(): 
                That set the layout of the plt figures.
                Note that after the end of each plot it is reset to the
                default matplotlib prop.
                
            loadData() :
                Import ENERGIES data  in a variable of type np.array.
                Note: data should be load before any attempting any
                operation on it. By default the the np.loadtxt is used for
                that. But there is an option to turn off the defaults and
                load data with more controls.
                
            repairTime() :
                Fix time problem when forget to activate accumulators in
                CPMD simulations. 
                
            plotNVE() :
                Plot NVE, NVT or NPT output ENERGIES.
                Note: for NPT, just have to turn on NPT in argumen, i.e
                NPT=True.
                
            plotANN() :
                Plot Annealing output ENERGIES.
                
            svfg() :
                To save figure when the automatic save is off.
                
        """
        print(msg)
        return print("Hope this doc"+
                     " :book: helped you :smiling_face_with_smiling_eyes: !"
                     )
        
    def changePath(self, newPath):
        """
        change the path to the complete path to a new working directory

        Parameters
        ----------
        newPath : str
            the new path.

        Returns
        -------
        str
            Print message indicating the new path.

        """
        self.path = newPath
        chdir(self.path)
        return print("path has changed to [red reverse]%s[/]" %self.path)
    
    @setDocStrAndApply  
    def setMatPltLibDefaults(self, default=False, style="seaborn-whitegrid",
                             rc_fig=dict(autolayout=True),
                             rc_ax=dict(labelweight="bold",
                                        labelsize="large",
                                        titleweight="bold",
                                        titlesize=14,
                                        titlepad=15,),
                             custom_cycler=(cycler(color=['c', 'm','y', 'k'])+
                                            cycler(lw=[1, 2, 3, 4])),
                             font = {'family' : 'sans-serif',
                                     'weight' : 'bold'}
                              ):

        matplotlib_params(default=default, style=style, rc_ax=rc_ax,
                          rc_fig=rc_fig, custom_cycler=custom_cycler,
                          font=font
                          )
        return print("plt default parameters have been setted")

    
    def repairTime(self):
        """
        Fix time problem when forget to activate accumulators in CPMD
        simulations.
        
    
        Parameters
        ----------
        time : np.array
            Input raw time.
    
        Returns
        -------
        np.array
            Time problem fixed.
    
        """
        self.data["time"] = repair_time(self.data["time"])
    
    def loadData(self, ncols=6, filename='ENERGIES', com="#", ignore=None,
                 sep=' ', skip=None, repTime=False, defaultLoad=True,
                 strcutDat=True, keepRawdat=False, NPT=False, meanT=False):
        """
         ata  in a variable of type array

        Parameters
        ----------
        filename : string
            Name of the file.
            The default is 'ENERGIES'
        ncols : int
            The total number of columns in data file.
            The default is 6.
        com : string, optional
            The comments symbol. Indicate line to be ignired.
            The default is "#".
        ignore : int, optional
            The number of line to be ignored. The default is None.
        sep : string, optional
            Columns seprator character.
            The default is ' '.
        skip : int, list or tuple, optional
            If specified the column number(s)  skip will be ignored.
            for e.g.: skip=0 then the first column will be ignored.
                      skip=[0, 1] or (0, 1) the two first columns we be
                      ignored.           
            The default is None.
        repTime : bool, optional
            If time (accumulator problem during MD simulation) to be repaired.
            The default is False.    
        defaultLoad : bool, optional
            load the data by default (with np.loadtxt). Else you have to load
            it after with the method loadData.
            The default is True

        Returns
        -------
        np.array
            array of the data imported.

        """
        if defaultLoad:
            self.data = np.loadtxt('ENERGIES', unpack=True)
        else:
            self.data = importData(filename, ncols=ncols, com="#",
                                       ignore=None, sep=' ', skip=None,
                                       repTime=False)
         
        if keepRawdat:
            self.raw_data=self.data.copy()
        
        self.data = structInpData2(self.data, NPT=NPT)
        if meanT:
            self.data['meanT']=_meanT(self.data.Tions)
        
    
    def plotNVE(self, savefig=True, NPT=False, insertEfic=True,
                posEfic=[0.30, 0.16, 0.10, .10],
                fig_param=dict(figsize=(16.4, 10.5), layout='constrained',
                               edgecolor='white', facecolor='white', dpi=160,
                               sharex=True
                               ),
                plot_param =dict(linewidth=3.5),
                plotmTemp=False, setAxLim:tuple=None
                ):
        """
        Plot for NVE, NVT or NPT output ENERGIES        

        Parameters
        ----------
        savefig : bool, optional
            If true the figure is directly saved. The default is True.
        NPT : bool, optional
            Indicate if NPT calculation. The default is False.
        insertEfic : bool, optional
            If True a supplementary ax is added. On this ax the fictitious
            energy is plotted in Hartree unit.
            The default is True.
        posEfic : list, optional
            The two first values give the position of the lower left corner.
            The two last are respectively the horizontal and the vertical sizes
            The default is [0.30, 0.16, 0.10, .10].
        fig_param : dict, optional
            parameters controlling the features of the figure, i.e the layout.
            The default is dict(figsize=(16.4, 10.5), layout='constrained',
                            edgecolor='white', facecolor='white', dpi=160 ).
        plot_param : dict, optional
            parameters controling the features of plot.
            The default is dict(linewidth=3.5).
        setAxLim : tuple, optional
            Set ylim of the ax that names is given as first element of the
            tuple and ylim being given as the remaining two element of the
            tuple.
            Possible ax names are :
                'Cons', "Vari"
                "Fict", "Temp"
                
            The default is None.
        Returns
        -------
        None.

        """
        time = self.data["time"] * self.tstep
        ax_names = ['Cons', 'Vari', 'Fict', 'Temp']
        plot_lab, plot_comp = setPlotLabels(NPT=NPT)
        
        #self.data["time"] = self.data["time"] * self.tstep
        self.data["DeltaEcons"] = (self.data.Econs - 
                                   self.data.Econs[0]) / self.data.Econs[0] 
        self.data["EkinIons"] = self.data.Tions * conv_kin(self.NbIons)
        self.data["Tfic"] = self.data.EkinFic  * conv_fic(self.NbElec)
        
        # settings for the plot
        _ax = {}
        fig, ax = plt.subplots(2, 2, **fig_param)
        
        _ax['Cons']=ax[0, 0]; _ax['Vari']=ax[0, 1]; _ax['Fict']=ax[1, 0]
        _ax['Temp']=ax[1, 1]
        # parameters for the legend of the different axis
        dic_legend = {}
        dic_legend['Cons'] = dict(loc=0, shadow=False, fontsize=20)
        dic_legend["Vari"] = dict(loc=0, shadow=False, fontsize=20)        
        dic_legend["Fict"] = dict(loc='center left', shadow=False, fontsize=20)        
        dic_legend["Temp"] = dict(loc=0, shadow=False, fontsize=25)

        #parameter for the labels of the different axis
        ax_lab = {} 
        ax_lab["Cons"] = (u'', u'(hartree)')
        ax_lab["Fict"] = (u'time (ps)', u'(hartree)')
        ax_lab["Vari"], ax_lab["Temp"] = (' ',' '), (u'time (ps)', u'(K)')
        
        for name in ax_names:
            for col in plot_comp[name]:
                # plotting of the axis
                _ax[name].plot(time, self.data[col],
                               label=plot_lab[col], **plot_param)
                # setting of the label to the axis
                _ax[name].set_xlabel(ax_lab[name][0], fontsize=20)
                _ax[name].set_ylabel(ax_lab[name][1], fontsize=20)
                #xticks
                _ax[name].tick_params(labelcolor='black', labelsize=18,
                                      width=3)
                # setting of the legend
                _ax[name].legend(**dic_legend[name])
                
        # insert e_fic plot inot the  ax[1,0]
        if insertEfic:
            inser = fig.add_axes(posEfic, facecolor='white')
            inser.plot(time, self.data.EkinFic)
            inser.set(title=u'fictious energy', xlabel=u'time (ps)',
                      ylabel='(hartree)')
        
        if plotmTemp:
            # Adding of additional informations on the plots
            _ax["Temp"].plot(time, self.data.meanT,
                              label=r'smoothed T$_ions$', alpha=0.6)
    
        
        # other operations
        plt.autoscale(True,'both')
        _ax["Vari"].set_ylim(-5e-4,5e-4)
        
        if setAxLim != None:
            _ax[setAxLim[0]].set_ylim(*setAxLim[1:])
            
        if savefig:
            plt.savefig('Energies.png')
            plt.show()
            plt.close(fig)
            
        else:
            plt.show()

            plt.close(fig)
            
    def plotANN(self, svfig=True, fig_param=dict(figsize=(12., 9.0),
                                                 layout='constrained',
                                                 edgecolor='white',
                                                 facecolor='white', dpi=160,
                                                 linewidth=0.0
                                                 ),
                plot_param=dict(linewidth=3.0)
                ):

        """
        Plot Annealing output ENERGIES        

        Parameters
        ----------
        svfig : bool, optional
            If true the figure is directly saved. The default is True.
        fig_param : dict, optional
            parameters controlling the features of the figure, i.e the layout.
            The default is dict(figsize=(12., 9.0), layout='constrained',
                            edgecolor='white', facecolor='white', dpi=160 ).
        plot_param : dict, optional
            parameters controling the features of plot.
            The default is dict(linewidth=3.0).

        Returns
        -------
        None.

        """        
        time = self.data["time"] * self.tstep # ps
        
        localisFig_tab=[['up','up'],['lowleft', 'lowright']]
        
        # creation of the figure object
        fig, ax = plt.subplot_mosaic(localisFig_tab, **fig_param)
        
        # data
        data_name ={}
        data_name["up"] = ('Econs', 'Eclas', "Epot")
        data_name["lowleft"] = ("EkinFic",)
        data_name["lowright"] = ("Tions",)
        
        # parameter for the plot
        #labels of the axis
        dic_plot ={}
        dic_plot["up"] = (r"$E_{cons}$", r"$E_{clas}$", r"$E_{pot}$")
        dic_plot["lowleft"] = (r"$E_{fic}$",)
        dic_plot["lowright"] = (r"$T_{ions}$",)
        
        # parameters for the legend
        dic_legend = {}
        dic_legend["up"] = dict(loc='upper right', shadow=True, fontsize=20)
        dic_legend["lowleft"] = dict(loc='upper right', shadow=True,
                                     fontsize=20)        
        dic_legend["lowright"] = dict(loc='upper right', shadow=True,
                                      fontsize=20)     
    
        #parameter for axis labels
        ax_lab = {} 
        ax_lab["up"]= ax_lab["lowleft"] =  (u'time (ps)', u'(hartree)')
        ax_lab["lowright"] = (u'time (ps)', u'(K)')
        
        #plotting
        for key in dic_plot.keys():
            for obs, lab in zip(data_name[key], dic_plot[key]):  
                ax[key].plot(time, self.data[obs], label=lab,
                             **plot_param)
                # setting of the label to the axis
                ax[key].set_xlabel(ax_lab[key][0], fontsize=20)
                ax[key].set_ylabel(ax_lab[key][1], fontsize=20)
                #xticks
                ax[key].tick_params(labelcolor='black', labelsize=18,
                                      width=3)
                # setting of the legend
                ax[key].legend(**dic_legend[key])
        
        # test
    
        # other operations
        if svfig:
            plt.autoscale(True,'both')
            plt.savefig('Energies.png')
            plt.show()
            plt.close(fig)
            
        else :
            plt.show()
            
            
    def svfg(self, fig):
        plt.autoscale(True,'both')
        plt.savefig('Energies.png')
        plt.show()
        plt.close(fig)
        
###############################################################################
#                                     AEMD                                    #
###############################################################################
class AEMD():
    """
        Class containing all the important plot in analysing AEMD simulation.
    """
    
    def __init__(self, path, NbIons=0, NbAlk=0, NbR=0, L=1., cellVol=1.,
                 alkr=False, inplane=True, phase2=False):
        """
        

        Parameters
        ----------
        path : str
            path to the directory of interest.
        NbIons : int, optional
            Total number of atoms in the system. The default is 0.
        NbAlk : int, optional
            Total number of atoms in alkyl chains. The default is 0.
        cellVol : float, optional
            The volume of the simulation box (in Angstrom cube).
            The default is 1.
        L : float, optional
            Spatially periodicity in the heat flux direction (in Angstrom).
            The default is 1.
        NbR : int, optional
            Total number of atoms in rings. The default is 0.
        
        alkr : bool, optional
            Precise wether or not alkyl subsystem is considered individually.
            The default is False.
        inplane : bool, optional
            Precise if the heat flux is inplane or not.
            The default is False.
        phase2 : bool, optional
            indicate wether of not the AEMD data to be plotted is for phase 2.
            In that case no fitting will be performed.
            The default is False.

        Returns
        -------
        None.

        """
        self.path=path
        chdir(path)
        self.AlkR=alkr
        self.lambd, self.volCell = L*1e-10, cellVol*1e-30
        self.inplan, self.phase = inplane, phase2
        self.N_ions, self.N_Alk, self.N_R = NbIons, NbAlk, NbR
        print("[bold magenta][cyan]*[/][/]"*79,
              "[bold magenta][cyan]*[/][/]"+ ' '*36 +
              '[bold magenta]AEMD[/]'+
              ' '*37 + "[bold magenta][cyan]*[/][/]",
              "[bold magenta][cyan]*[/][/]"*79, sep='\n') 
        print(
            """
            For informations about all the available methods :
                -> Apply getInfos()
            For informations about the arguments of a given method :
                -> Apply infoMeth(obj.method); with obj the object you created
                    and its method.
                    
                    Or simply use :
                        >> print(instance.method.__doc__)

                            
            """
            )

    def infoMeth(self, method):
        """
        Print the documenation about a given method.

        Parameters
        ----------
        method : method
            The method of interest.
            e.g: if func() is a method of the class object h / h = Test()
            h.infoMeth(h.func) print the documentation of the method func()
        Returns
        -------
        None.

        """
        print(method.__doc__)
    
    def getInfos(self):
        msg ="""
        This Class Contains the following methods:
            changePath() :
                If you want to change the path from the one you
                defined when creating the object structuralProp().
            
            setMatPltLibDefaults(): 
                That set the layout of the plt figures.
                Note that after the end of each plot it is reset to the
                default matplotlib prop.
            
            repairTime() :
                Fix time problem when forget to activate accumulators in
                CPMD simulations.
                
            loadData() :
                Load data of temperature trajectory.
                Note: data should be load before trying any
                plot or operation. The only execption is for the method 
                plotblc2Alkblc1Ring()
                
            plot4totSystm() :
                Plot evolution of the temperatures (bloc1 and bloc2) of
                the total system.
                
            plot4Alkyls() :
                Plot temperature evolution toward equiibrium for alkyls
                chains
                
            plot4Rings() :
                Evolution of the Rings temperatures (hot bloc
                and cool bloc) towards equilibrium.
                
            profilT() :
                Plot the profil of the temperature along the heat flux
                direction.
                
            mergePhase() :
               merge phase1 temperature trajectory with phase 2
                
            plotblc2Alkblc1Ring() :
                Plot the evolution of both the temperature of rings and alkyls.
                This plot is only relevant when interested in heat exchange
                between them.
        """
        print(msg)
        return print("Hope this doc"+
                     " :book: helped you :smiling_face_with_smiling_eyes: !"
                     )
        
    def changePath(self, newPath):
        """
        change the path to the complete path to a new working directory

        Parameters
        ----------
        newPath : str
            the new path.

        Returns
        -------
        str
            Print message indicating the new path.

        """
        self.path = newPath
        chdir(self.path)
        return print(":warning:path has changed to"+
                     " [red reverse]%s[/]" %self.path)
    
    
    @setDocStrAndApply  
    def setMatPltLibDefaults(self, default=False, style="seaborn-whitegrid",
                             rc_fig=dict(autolayout=True),
                             rc_ax=dict(labelweight="bold",
                                        labelsize="large",
                                        titleweight="bold",
                                        titlesize=14,
                                        titlepad=15,),
                             custom_cycler=(cycler(color=['c', 'm','y', 'k'])+
                                            cycler(lw=[1, 2, 3, 4])),
                             font = {'family' : 'sans-serif',
                                     'weight' : 'bold'}
                              ):

        matplotlib_params(default=default, style=style, rc_ax=rc_ax,
                          rc_fig=rc_fig, custom_cycler=custom_cycler,
                          font=font
                          )
        return print("plt default parameters have been setted")
        
            
    def loadData(self, path=".\\", inplace=True, tstep:int=None,
                 reset_orgin_time=False):
        """
        Load data of temperature trajectory.
        path : str, optional
            path to the data
            The default is the cwd.
        inplace: bool, optional
            if True, the data are autically attributed. Else if return data.
            The default is True
        Returns
        -------
        if AlkR and not inplace :
            return time, T_bloc1, T_bloc2, time2, T_Alk1, T_Alk2, T_R1, T_R2
        elif AlkR and not inplace:
            time, T_bloc1, T_bloc2

        """
        if self.AlkR:
            if inplace:
                self.time, self.T_bloc1, self.T_bloc2 = np.loadtxt(
                                                            path+'trajDT.dat',
                                                            unpack=True
                                                            )
                
                self.time2, self.T_Alk1, self.T_Alk2 = np.loadtxt(
                                                        path+'trajDTAlk.dat',
                                                        unpack=True
                                                        )
                _, self.T_R1, self.T_R2 = np.loadtxt(
                                                    path+'trajDTRings.dat',
                                                    unpack=True)
            else:           
                time, T_bloc1, T_bloc2 = np.loadtxt(path+'trajDT.dat',
                                                    unpack=True)
                time2, T_Alk1, T_Alk2 = np.loadtxt(path+'trajDTAlk.dat',
                                                   unpack=True)
                _, T_R1, T_R2 = np.loadtxt(path+'trajDTRings.dat',
                                           unpack=True)
                return time,T_bloc1, T_bloc2, time2, T_Alk1, T_Alk2, T_R1, T_R2
            if reset_orgin_time:
                self.time2 = self.time2 - self.time2[0]
            if tstep != None:
                self.time2 = self.time2*tstep*0.024*1e-3
                
                
        else:
            if inplace:
                self.time, self.T_bloc1, self.T_bloc2 = np.loadtxt(
                                                            path+'trajDT.dat',
                                                               unpack=True)
            else:
                time, T_bloc1, T_bloc2 = np.loadtxt(path+'trajDT.dat',
                                                    unpack=True
                                                    )
                return time, T_bloc1, T_bloc2
        if reset_orgin_time:
            self.time = self.time - self.time[0]
        if tstep != None:
            self.time = self.time*tstep*0.024*1e-3

    def repairTime(self):
        """
        Fix time problem when forget to activate accumulators in CPMD
        simulations.
        
    
        Parameters
        ----------
        time : np.array
            Input raw time.
    
        Returns
        -------
        np.array
            Time problem fixed.
    
        """
        if self.AlkR:
            self.time = repair_time(self.time)
            self.time2 = repair_time(self.time2)
        else:
            self.time = repair_time(self.time)        
        
    def mergePhase(self, path2, tstep:int=None):
        """
        merge phase1 temperature trajectory with phase2
        

        Parameters
        ----------
        path2 : str
            Path to the temperature trajectory in phase 1.

        Returns
        -------
        None.

        """
        if self.AlkR:
            data = self.loadData(path=path2, inplace=False, tstep=tstep)
            time, T_bloc1, T_bloc2, time2, T_Alk1, T_Alk2, T_R1, T_R2 = data
            if tstep != None:
                time = time * tstep * 0.024*1e-3
                time2 = time2 * tstep * 0.024*1e-3
            self.Mergedtime = np.concatenate((time, self.time + time[-1]))
            self.Mergedtime2 = np.concatenate((time, self.time2 +time[-1]))
            self.MergedT_bloc1 = np.concatenate((T_bloc1, self.T_bloc1))
            self.MergedT_bloc2 = np.concatenate((T_bloc2, self.T_bloc2))
            self.MergedT_Alk1 = np.concatenate((T_Alk1, self.T_Alk1))
            self.MergedT_Alk2 = np.concatenate((T_Alk2, self.T_Alk2))
            self.MergedT_R1 = np.concatenate((T_R1, self.T_R1))
            self.MergedT_R2 = np.concatenate((T_R2, self.T_R2))
                
        else:
            time, T_bloc1, T_bloc2 = self.loadData(path=path2, inplace=False)
            if tstep != None:
                time = time * tstep * 0.024*1e-3
            self.Mergedtime = np.concatenate((time, self.time + time[-1]))
            self.Mergedtime2 = np.concatenate((time, self.time2 +time[-1]))
            self.MergedT_bloc1 = np.concatenate((T_bloc1, self.T_bloc1))
            self.MergedT_bloc2 = np.concatenate((T_bloc2, self.T_bloc2))
        
        
                    
    def plot4totSystm(self, mergedData=False, onlyTemp=True,figs=(9, 7.5),
                      timeLim=False,
                      endTime=None, delay=False, wait=None, setlimDT=False,
                      limDT=(10,300), xlim=None, posText=(.01, 0.1),
                      separate=False
                      ):
        """
        plot evolution of the temperatures (bloc1 and bloc2) of the total
        system.

        Parameters
        ----------
        mergedData : bool, optional,
            Indicate if data you want to see if a merged data (of phase 1
                                                               and 2).
            The default is False.
            
            Note : if True, the merging must be done before by using
                    mergePhase() method
        timeLim : bool, optional
            Activate the option to fix a limit to the time for plotting and
            other analysis. The default is False.
        endTime : float, optional
            this time limit (must be fixed as float if TtimeLim is True).
            (filtre passe haut)
            The default is None.
        delay : bool, optional
            Activate the filter that consist in starting the plot latter.
            The default is False.
            (filtre passe bas)
        wait : float, optional
            All data such as the time is inferior to wait are not take into
            account. The default is None.
        setlimDT : bool, optional
            Indicate if to set limits of the temperature difference
            as function of the time.
            The default is False.
            (filtre passe bas)
        limDT : tuple, optional
            tuple of the limit of the temperature difference.
            The default is (10, 300)
        xlim : float, optional
            val of xlim (only the plot is zoomed, data is not affected).
                         The default is None.
        posText : tuple, optional
            Position of the text.
            The first value is the position in x, the second in y and the last
            for interline distance.
            The default is (.01, 0.1).

        Returns
        -------
        None.

        """

        if not mergedData:
            time, Tb1, Tb2 = opTemp(self.time, self.T_bloc1, self.T_bloc2,
                   timeLim=timeLim, endTime=endTime, delay=delay, wait=wait)
         
                # plot AEMD p3ht
            fig, ax = plot_AEMD(time, Tb1, Tb2, 'trajDT.png', title='P3HT',
                                posi=posText, N=self.N_ions, xlim=xlim,
                                inplane=self.inplan, Lm=self.lambd,
                                cellVol=self.volCell,
                                setlimDT=setlimDT, limDT=limDT,
                                phase2=self.phase, separate=separate,
                                svnameDeltaT="DeltaT.png"
                                )
            #plt.rcdefaults()
            return fig, ax
        else:
            time, Tb1, Tb2 = opTemp(self.Mergedtime, self.MergedT_bloc1,
                                    self.MergedT_bloc2,
                   timeLim=timeLim, endTime=endTime, delay=delay, wait=wait)
         
                # plot AEMD p3ht
            if onlyTemp:
                fig, ax = plt.subplots(figsize=figs, dpi=300)
                ax.plot(time, Tb1, label=r'T$_\mathrm{C}$', color='tomato')
                ax.plot(time, Tb2, label=r'T$_\mathrm{F}$', color='c')
                ax.legend(loc=0)
                ax.set_xlabel(r"Time (ps)")
                ax.set_ylabel(r'T (K)')
                
                if xlim is not None:
                    ax.set_xlim(xlim)
                    ax.set_xlim(xlim)
                else:
                    plt.autoscale(True,'both')
                if setlimDT:
                    ax.set_ylim(limDT)
                    
                # other operations
                plt.savefig("mergedTrajDT.png", dpi=300)
                plt.show()
                plt.close(fig)
            else:
                fig, ax = plot_AEMD(time, Tb1, Tb2, 'trajDT.png', title='P3HT',
                                    posi=posText, N=self.N_ions, xlim=xlim,
                                    inplane=self.inplan, Lm=self.lambd,
                                    cellVol=self.volCell,
                                    setlimDT=setlimDT, limDT=limDT,
                                    phase2=False
                                    )
            #plt.rcdefaults()
            return fig, ax
        
    def plot4Alkyls(self, mergedData=False, timeLim=False, endTime=None,
                    delay=False, wait=None,
              setlimDT=False, limDT=(10,300),
              xlim=None, posText=(.01, .10), 
              separate=False
              ):
        """
        plot temperature evolution toward equiibrium for alkyls chains

        Parameters
        ----------
        mergedData : bool, optional,
            Indicate if data you want to see if a merged data (of phase 1
                                                               and 2).
            The default is False.
            
            Note : if True, the merging must be done before by using
                    mergePhase() method
        timeLim : bool, optional
            Activate the option to fix a limit to the time for plotting and
            other analysis. The default is False.
        endTime : float, optional
            this time limit (must be fixed as float if TtimeLim is True).
            (filtre passe haut)
            The default is None.
        delay : bool, optional
            Activate the filter that consist in starting the plot latter.
            The default is False.
            (filtre passe bas)
        wait : float, optional
            All data such as the time is inferior to wait are not take into
            account. The default is None.
        setlimDT : bool, optional
            Indicate if to set limits of the temperature difference
            as function of the time.
            The default is False.
            (filtre passe bas)
        limDT : tuple, optional
            tuple of the limit of the temperature difference.
            The default is (10, 300).
        xlim : float, optional
            val of xlim (only the plot is zoomed, data is not affected).
                         The default is None.
        posText : tuple, optional
            Position of the text.
            The first value is the position in x, the second in y.
            The default is (.01, .10).

        Returns
        -------
        None.

        """
        if not mergedData:
            time, Tb1, Tb2 = opTemp(self.time, self.T_Alk1, self.T_Alk2,
                   timeLim=timeLim, endTime=endTime, delay=delay, wait=wait)
         
            # plot AEMD p3ht
            fig, ax = plot_AEMD(time, Tb1, Tb2, 'trajDTAlk.png',
                                        title='Alkyls', posi=posText,
                                        N=self.N_Alk, xlim=xlim,
                                        inplane=self.inplan, Lm=self.lambd,
                                        cellVol=self.volCell,
                                        setlimDT=setlimDT, limDT=limDT,
                                        phase2=self.phase, separate=separate,
                                        svnameDeltaT="DeltaTAlk.png"
                                        )
            return fig, ax
        else:
            time, Tb1, Tb2 = opTemp(self.Mergedtime2, self.MergedT_Alk1,
                                    self.MergedT_Alk2,
                   timeLim=timeLim, endTime=endTime, delay=delay, wait=wait)
         
                # plot AEMD p3ht
            fig, ax = plot_AEMD(time, Tb1, Tb2, 'trajDTAlk.png',
                                        title='Alkyls', posi=posText,
                                        N=self.N_Alk, xlim=xlim,
                                        inplane=self.inplan, Lm=self.lambd,
                                        cellVol=self.volCell,
                                        setlimDT=setlimDT, limDT=limDT,
                                        phase2=False
                                        )
            return fig, ax
        
    def plot4Rings(self, mergedData=False, timeLim=False, endTime=None,
                   delay=False, wait=None,
                   setlimDT=False, limDT=(10,300),
                   xlim=None, posText=(.01, .10), 
                   separate=False
                   ):
        """
        Evolution of the Rings temperatures (hot bloc and cool bloc) towards
        equilibrium. 

        Parameters
        ----------
        mergedData : bool, optional,
            Indicate if data you want to see if a merged data (of phase 1
                                                               and 2).
            The default is False.
            
            Note : if True, the merging must be done before by using
                    mergePhase() method
        timeLim : bool, optional
            Activate the option to fix a limit to the time for plotting and
            other analysis. The default is False.
        endTime : float, optional
            this time limit (must be fixed as float if TtimeLim is True).
            (filtre passe haut)
            The default is None.
        delay : bool, optional
            Activate the filter that consist in starting the plot latter.
            The default is False.
            (filtre passe bas)
        wait : float, optional
            All data such as the time is inferior to wait are not take into
            account. The default is None.
        setlimDT : bool, optional
            Indicate if to set limits of the temperature difference
            as function of the time.
            The default is False.
            (filtre passe bas)
        limDT : tuple, optional
            tuple of the limit of the temperature difference.
            The default is (10, 300).
        xlim : float, optional
            val of xlim (only the plot is zoomed, data is not affected).
                         The default is None.
        posText : tuple, optional
            Position of the text.
            The first value is the position in x, the second in y and the last
            for interline distance.
            The default is (.01, .10).

        Returns
        -------
        None.

        """
        if not mergedData:
            time, Tb1, Tb2 = opTemp(self.time, self.T_R1, self.T_R2,
                   timeLim=timeLim, endTime=endTime, delay=delay, wait=wait)
         
                # plot AEMD p3ht
            fig, ax = plot_AEMD(time, Tb1, Tb2, 'trajDTRings.png',
                                        title='Rings', posi=posText,
                                        N=self.N_R, xlim=xlim,
                                        inplane=self.inplan,
                                        Lm=self.lambd,
                                        cellVol=self.volCell,
                                        setlimDT=setlimDT, limDT=limDT,
                                        phase2=self.phase, separate=separate,
                                        svnameDeltaT="DeltaTRings.png"
                                        )
        else:
            time, Tb1, Tb2 = opTemp(self.Mergedtime2, self.MergedT_R1,
                                    self.MergedT_R2,
                                    timeLim=timeLim, endTime=endTime,
                                    delay=delay, wait=wait)
         
                # plot AEMD p3ht
            fig, ax = plot_AEMD(time, Tb1, Tb2, 'trajDTRings.png',
                                        title='Rings', posi=posText,
                                        N=self.N_R, xlim=xlim,
                                        inplane=self.inplan,
                                        Lm=self.lambd,
                                        cellVol=self.volCell,
                                        setlimDT=setlimDT, limDT=limDT,
                                        phase2=False
                                        )
            
    def plotblc2Alkblc1Ring(self, path=".\\", timeLim=False, endTime=None,
                            delay=False, wait=None, setlimDT=False,
                            limDT=(10,300), xlim=None, posText=(.01, .10)
                            ):
        """
        Plot the evolution of both the temperature of rings and alkyls.
        This plot is only relevant when interested in heat exchange between
        the two.

        Parameters
        ----------
        path : str, optional
            path to the data
            The default is the cwd.
        timeLim : bool, optional
            Activate the option to fix a limit to the time for plotting and
            other analysis. The default is False.
        endTime : float, optional
            this time limit (must be fixed as float if TtimeLim is True).
            (filtre passe haut)
            The default is None.
        delay : bool, optional
            Activate the filter that consist in starting the plot latter.
            The default is False.
            (filtre passe bas)
        wait : float, optional
            All data such as the time is inferior to wait are not take into
            account. The default is None.
        setlimDT : bool, optional
            Indicate if to set limits of the temperature difference
            as function of the time.
            The default is False.
            (filtre passe bas)
        limDT : tuple, optional
            tuple of the limit of the temperature difference.
            The default is (10, 300).
        xlim : float, optional
            val of xlim (only the plot is zoomed, data is not affected).
                         The default is None.
        posText : tuple, optional
            Position of the text.
            The first value is the position in x, the second in y and the last
            for interline distance.
            The default is (.01, .10).

        Returns
        -------
        None.

        """
        time, T_R, T_A = np.loadtxt("trajDTAlRi.dat", unpack=True)
        
         
                # plot AEMD p3ht
        fig, ax = plot_AEMD(time, T_R, T_A, 'trajDTRings.png',
                                    title='', posi=posText,
                                    N=self.N_R, xlim=xlim,
                                    inplane=self.inplan,
                                    Lm=self.lambd,
                                    cellVol=self.volCell,
                                    setlimDT=setlimDT, limDT=limDT,
                                    phase2=False,
                                    titlBlc1="Ring", titleBlc2="Alkyls"
                                    )

    def profilT(self,
            posText=(.01, .10), h_alignment='left',figsize=(6.4, 4.8),
            fig_param = dict(layout='constrained',
                             edgecolor='white', facecolor='white',
                             dpi=150, linewidth=2.0),
            plot_param = dict(marker=None, linestyle=None,
                              linewidth=3.0, markersize=None),
            translate_slice=False, min_val_slice:float=0., max_val_slice=0.,
            path=None, axis=None, figr=None, svname="profilT.png",
            label='Phase 1', color:str=None, color_fit:str=None,
            reverse:bool=False
            ):
        """
        Plot the profil of the temperature along the heat flux direction.

        Parameters
        ----------
        posText : tuple, optional
            Position of the text that will be added is phase 2.
            first value is the position in x.
            the second is the positiion in y.
            The default is (.01, .1).
        h_alignment : str, optional
            Horizontal alignment of the texts
            The default is "left"
        figsize: tuple, optional
            The default is (6.4, 4.8).
        fig_param : dict, optional
            Contains parameters controlling the plot layout characteristics.
            The default is dict(layout='constrained',
                                edgecolor='white', facecolor='white',
                                dpi=150, linewidth=2.0).
        plot_param : dict, optional
            Contains parameters controlling the characteristic of the plot.
            The default is dict(marker=None, linestyle=None,
                                linewidth=3.0, markersize=None).

        Returns
        -------
        fig : plt.figure
        ax : plt.axis

        """
        if path !=None:
            _path=path
        else:
            path=''
        try:
            open(path+'profilT.data')
        except:
            dist, temp, *_ = np.loadtxt(path+'profilT.dat', unpack=True)
        else:
            dist, temp, *_ = np.loadtxt(path+'profilT.data', unpack=True)
        if translate_slice:
            temp=reorganize_data(dist, temp, min_val=min_val_slice,
                                   maxval=max_val_slice
                                   )

        x_=dist
        y_=temp
        if axis == None and figr == None:
            #reorganize_data(dist, 2)
            # settings for the plot
            if reverse:
                figsize = figsize[-1], figsize[0]
                
            fig, ax = plt.subplots(nrows=1, ncols=1,figsize=figsize,
                                   **fig_param)

            #parameter for the labels of the different axis
            if self.inplan:
                ax_lab = (r'x ($\AA$)', r'Température (K)')
            else:
                ax_lab = (r'z ($\AA$)', r'Température (K)')
            # setting of the label to the axis

            #ticks
            ax.tick_params(labelcolor='black', labelsize=18, width=3)
            x_lab = ax_lab[0]
            y_lab = ax_lab[1]
            if reverse:
                x_, y_ = y_, x_
                x_lab, y_lab = y_lab, x_lab
                
            ax.set_xlabel(x_lab, fontsize=20)
            ax.set_ylabel(y_lab, fontsize=20)

            #plotting 
            ax.plot(x_, y_, **plot_param,
                    label=label, color=color)
        
            if self.phase and self.inplan:
                # fit par func sin
                funcSin = lambda x, A, B, phi: A*np.sin(
                                                    2*(np.pi/(self.lambd*1e10))*x +
                                                            phi) + B
                popt, pcov = curve_fit(funcSin, dist, temp)
                x_fit = dist
                y_fit = funcSin(dist,*popt)
                if reverse:
                    x_fit, y_fit = y_fit, x_fit
                    
                ax.plot(x_fit, y_fit,
                        label=r'fit $A \times \sin{(\frac{2\pi}{L_x} x + \phi)} + B$',
                        color=color_fit, lw=2.5)
                
                # Adding of additional informations on the plots
                ax.text(*posText, r'A = {:.2f} K ; B = {:.2f} K'.format(
                    popt[0], popt[1]), fontsize=10, color='k',
                    transform=ax.transAxes, horizontalalignment=h_alignment)
                ax.text(posText[0], posText[1]+.10, r'$\phi$ = %.2f' %popt[2],
                        fontsize=10, color='k', transform=ax.transAxes,
                        horizontalalignment=h_alignment)
            
    
            
    
            
            # other operations
            plt.autoscale(True,'both')
            plt.savefig(svname)
            plt.show()
            plt.close(fig)
        else:
            # parameters for the legend of the different axis
            dic_legend= dict(loc=0, shadow=False, fontsize=15)
            if reverse:
                x_, y_ = y_, x_
                
            axis.plot(x_, y_, **plot_param,
                    label=label, color=color)
            axis.legend(**dic_legend)
            figr.savefig(svname)
            fig = figr
            ax = axis
        return fig, ax

###############################################################################
#                                     DOS                                     #
###############################################################################
class DensityOfState:
    
    """
        Contains methods to plot the density of state
    """
    
    def __init__(self, filename, path, typeDos, cpmd=False):
        """
        

        Parameters
        ----------
        filename : str
            The name of the input file.
            Either the cpmd output or a file with eigen value and gaussian
            broadened val.
        path : str
            the complete path to the directory containing the file.
        typeDos : str, optional
            Only two possible values: 'e' or 'v'.
            'e' for electronic
            'v' for vibrational.
            The default is 'v'.
        cpmd : bool, optional
            if True, that indicate that eigenvalues will be extract directly
            from a cpmd output file. Otherwise the input file must exactly have
            two columns. The first corresponding to intensities and the last
            for the position in energy.
        """
        print("[bold magenta][cyan]*[/][/]"*79,
              "[bold magenta][cyan]*[/][/]"+' '*30+
              "[bold magenta]DENSITY OF STATES[/]"+
              ' '*30+"[bold magenta][cyan]*[/][/]",
              "[bold magenta][cyan]*[/][/]"*79,
              sep='\n') 
        self.path = path
        chdir(self.path)
        self.file = filename
        self.typeDos =typeDos
        if not cpmd:
            self.eigen, self.au = np.loadtxt(filename, unpack=True)
            if self.typeDos=='e':
                self.fermi = float(
                    input("Entrez la valeur en eV du niveau de fermi : ")
                    )
        print(
            """
            For informations about all the available methods :
                -> Apply getInfos()
            For informations about the arguments of a given method :
                -> Apply infoMeth(obj.method); with obj the object you created
                    and its method.
                                        
                    Or simply use :
                        >> print(instance.method.__doc__)
            """
            )
        
    def infoMeth(self, method):
        """
        Print the documenation about a given method.

        Parameters
        ----------
        method : method
            The method of interest.
            e.g: if func() is a method of the class object h / h = Test()
            h.infoMeth(h.func) print the documentation of the method func()
        Returns
        -------
        None.

        """
        print(method.__doc__)
        
    def getInfos(self):
        msg ="""
        This Class Contains the following methods:
            changePath() :
                If you want to change the path from the one you
                defined when creating the object structuralProp().
            
            setMatPltLibDefaults(): 
                That set the layout of the plt figures.
                Note that after the end of each plot it is reset to the
                default matplotlib prop.
                
            xtracEigenValFromCpmd() :
                Extract eigenvalues directly from cpmd output file.
                
            broadening() :
                guassian broadening of the eigenvalues spectrum.
                Note : That method should be applied before trying to plot
                the DOS.
                
            plotDOS() :
                Plot the density of the state (DOS).
                
        Note :
            There is no method loadData since data are automatically
            loaded.
                
        """
        print(msg)
        return print("Hope this doc"+
                     " :book: helped you :smiling_face_with_smiling_eyes: !"
                     )
        
    def changePath(self, newPath):
        """
        change the path to the complete path to a new working directory

        Parameters
        ----------
        newPath : str
            the new path.

        Returns
        -------
        str
            Print message indicating the new path.

        """
        self.path = newPath
        chdir(self.path)
        return print("path has changed to [red reverse]%s[/]" %self.path)

    @setDocStrAndApply  
    def setMatPltLibDefaults(self, default=False, style="seaborn-whitegrid",
                             rc_fig=dict(autolayout=True),
                             rc_ax=dict(labelweight="bold",
                                        labelsize="large",
                                        titleweight="bold",
                                        titlesize=14,
                                        titlepad=15,),
                             custom_cycler=(cycler(color=['c', 'm','y', 'k'])+
                                            cycler(lw=[1, 2, 3, 4])),
                             font = {'family' : 'sans-serif',
                                     'weight' : 'bold'}
                              ):

        matplotlib_params(default=default, style=style, rc_ax=rc_ax,
                          rc_fig=rc_fig, custom_cycler=custom_cycler,
                          font=font
                          )
        return print("plt default parameters have been setted")
        
    

    def xtracEigenValFromCpmd(self, NbOfUnOc:int=10):
        """
        Extract eigenvalues directly from cpmd output file.

        Parameters
        ----------
        NbOfUnOc : int, optional
            The number of unoccupied states to be taken into account.
            The default is 10.

        Returns
        -------
        None.

        """

        prekeyWord = dict(e='FINAL RESULTS',
                          v='PURIFICATION OF DYNAMICAL MATRIX'
                          )
        keyWords = dict(e='EIGENVALUES(EV) AND OCCUPATION:',
                        v='HARMONIC FREQUENCIES [cm**-1]:'
                        )
        endWords = dict(e='CHEMICAL POTENTIAL =',
                        v='ChkSum(FREQ) ='
                        )
        fileOut = dict(e='rawedos', v='rawvdos')
        self.fermi = dict(e=0, v=None)
        ncols = dict(e=6, v=4)
        skip = dict(e=False, v=True)
        
        self.data, self.fermi[self.typeDos] = extract(self.file,
                                                      fileOut[self.typeDos],
                                        prekeyWord[self.typeDos],
                                        keyword=keyWords[self.typeDos],
                                        endword=endWords[self.typeDos],
                                        skip=skip[self.typeDos],
                                        ncols=ncols[self.typeDos],
                                        fermi=self.fermi[self.typeDos],
                                        NbUnOcState=NbOfUnOc
                                        )

        
    def broadening(self, sigma=0.005, minVal=None, maxVal=None, Npoints=2000,
                   save:bool=False, savename='dos.txt'
                   ):
        """
        guassian broadening of the eigenvalues spectrum.

        Parameters
        ----------
        sigma : float, optional
            guassian width.
            The default is 0.005.
        minVal : float, optional.
            The minimum value x-axis. The default is None.
        maxval : float, optional.
            The maximum value. The default is None.
        Npoints : int, optional
            Total number of points at which to evaluate the gaussian
            contribution of the ensemble of the eigenvalues.
            The default is 2000.
        save : bool, optional
            If true, the result of the operation will be save as savename.
            The default is False.
        savename: str, optional
            The name of the file to be saved.

        """
        forma = '{:8.3f}    {:8.4f}\n'
        #occupation = []
        if minVal or maxVal is None:
            minVal, maxVal = min(self.data), max(self.data)
            
        self.eigen, self.au  = gaussBroadn(Eigen=self.data, minVal=minVal,
                                           maxVal=maxVal, npoints=Npoints,
                                           width=sigma
                                           )

        if save:
            with open(savename, "w") as f:
                for i, j in zip(self.eigen, self.au):
                    f.write(forma.format(i, j))
  
    def plotDOS(self, fig_param=dict(figsize=(6.4, 4.5),
                                           layout='constrained',
                                           edgecolor='white',
                                           facecolor='white',
                                           dpi=200,
                                           linewidth=0.0
                                           ),
                plot_param=dict(marker=None, linestyle='-',
                                linewidth=1, markersize=None
                                ),
                alpha=.5,
                colorFermi='red',
                limits=[(0, 3.65), (0, 500)],
                tick_param=dict(labelcolor='black', labelsize=18, width=3,
                                grid_linewidth=5, direction='out'
                                ),
                maxFrmiLine:float=5., noYticks=True,
                pathTosec:str=None
                ,shift2:float=0.0, fermi2:float=0.0,
                labels=('',), locLeg="upper left", bbox_leg=(1.05, 1.0)
                ):
        """
        Plot the density of state.
    
        Parameters
        ----------
        fig_param : dict, optional
            control the figure layout. The default is dict(figsize=(6.4, 4.5),
                                             layout='constrained',
                                             edgecolor='white',
                                             facecolor='white',
                                             dpi=160,
                                             linewidth=0.0 
                                             ).
        plot_param : dic, optional
            control the plot style.
            The default is dict(marker=None, linestyle='-',
                                linewidth=None, markersize=None).
        alpha : float, optional
            control the transparency of the filling. The default is .5.
        colorFermi : str, optional
            color of the fermi line plot. The default is 'red'.
        limits : list of tuple, optional
            the couple of xlim and ylim. The default is [(0, 3.65), (0, 500)].
        tick_param : dict, optional
            Control the axis ticks style.
            The default is tick_param=dict(labelcolor='black',
                                           labelsize=18, width=3)
        maxFrmiLine : float, optional.
            Theheight of the fermi line. The defaut is 5.0. 
        noYticks : bool, optional.
            if True then will not put yticks and label 
            The default is False.
        pathTosec: str, optional.
            Complete path to acces to the file of the second eDos to be plotted
            The default is None.
        fermi2 : float, optional.
            fermi level of the second dos    
            The default is 0.
        shift2 : float, optional
            Value of the vertical shift to be carried on the second doc.
            The default is 0.
        labels : tuple, optional
            Labels of the dos to be simultaneously plotted
            The default is ("",).
        locLeg : str, optional
            The location of the axes legend. 
            The default is "upper left"
        bbox_leg : tuple of float, optional
            The position of the legend texts. 
            The default is (1.05, 1.0).
        Returns
        -------
        None.
    
        """  
        
        # settin of figure and ax
        fig, ax = plt.subplots(**fig_param)
        
        # label for the different graph
        dic_plot = dict(e=r"e-DOS", v= r"vib-DOS")
        
        # parameters of the labels of the axes
        ax_lab = dict(e=r"Energies (eV)", v=r"Frequency $(cm^{-1})$")
  
        # plotting of the axis
        ax.plot(self.eigen, self.au,
                label=labels[0],
                **plot_param
                )
        
        
        # fill between plot
        ax.fill_between(self.eigen, self.au, where=self.au>=min(self.au),
                        alpha=alpha, linewidth=0
                        )
        
        #xticks
        ax.tick_params(**tick_param)
        
        fermi = self.fermi[self.typeDos]
        
        # Adding of additional informations on the plots
        if fermi is not None:
            ax.vlines(fermi, min(self.au), maxFrmiLine, colors=colorFermi,
                      linestyles='dashed', label='Niveau de Fermi'
                      )
        
        
        if pathTosec != None:
            eig, au = np.loadtxt(pathTosec, unpack=True)

            au += shift2
            ax.plot(eig, au,
                    label=labels[-1],
                    **plot_param
                    )

            ax.fill_between(eig, au, min(au),
                            where=au>=min(au),alpha=alpha, linewidth=0)
            ax.vlines(fermi2, min(au), maxFrmiLine+shift2,
                      colors=colorFermi,
                      linestyles='dashed')#, label='fermi level'
        if self.typeDos == 'e':
            ax.legend(bbox_to_anchor=bbox_leg, loc=locLeg)
            
            #ax.legend(loc=0)
        # setting of the label to the axis
        ax.set_xlabel(ax_lab[self.typeDos], fontsize=18)
        if noYticks:
            plt.yticks([])
        else:
            ax.set_ylabel(r"a. u.", fontsize=20)

    
        plt.autoscale("both")
        if isinstance(limits[0], tuple):
            ax.set_xlim(*limits[0])
        if isinstance(limits[1], tuple):
            ax.set_ylim(*limits[1])

        figname={}
        figname= dict(e='eDOS.png', v='vibDOS.png')
    
        plt.savefig(figname[self.typeDos], dpi=300)
        plt.show()
        plt.close(fig)

###############################################################################
#                                   STRUCTURE                                 #
###############################################################################
class StructuralProp:
    """
     Contains important plot concerning structural properties.
        Parameters
        ----------
        path : str
            path to the directory of interest..
    """
    def __init__(self, path):
        """
        

        Parameters
        ----------
        path : str
            path to the directory of interest.

        Returns
        -------
        None.

        """
        self.path = path
        chdir(self.path)
        
        print("[bold cyan]*[/]"*79,
              "[bold cyan]*[/]"+' '*28+
              '[bold magenta]STRUCTURAL PROPERTIES[/]'+
              ' '*28+"[bold cyan]*[/]",
              "[bold cyan]*[/]"*79, sep='\n')        
        print(
            """
            For informations about all the available methods :
                -> Apply getInfos()
            For informations about the arguments of a given method :
                -> Apply infoMeth(obj.method); with obj the object you created
                    and its method.
                                        
                    Or simply use :
                        >> print(instance.method.__doc__)
            """
            )
    
    def infoMeth(self, method):
        """
        Print the documenation about a given method.

        Parameters
        ----------
        method : method
            The method of interest.
            e.g: if func() is a method of the class object h / h = Test()
            h.infoMeth(h.func) print the documentation of the method func()
        Returns
        -------
        None.

        """
        print(method.__doc__)
    
    def getInfos(self):
        msg ="""
        This Class Contains the following methods:
            changePath() : If you want to change the path from the one you
                defined when creating the object structuralProp().
            
            setMatPltLibDefaults(): That set the layout of the plt figures.
                Note that after the end of each plot it is reset to the default
                matplotlib prop.
            
            plotBondAngleDist() :
                To plot the bond angle distribution
                
            plotPartialgOfR() :
                Plot the partial g(r).
                
            plotMSDAndDiffCoef() :
                Plot the MSD and the diffusion coefficient related to the
                structure.
                
            plotStress() :
                Plot the evolution of the stress.
                
            PlotStressAndCellDim() :
                Plot both the stres and the cell dimension as function
                the time.
        Note : there is no need to preload data, they will be loaded
                when needed.
        """
        print(msg)
        return print("Hope this doc"+
                     " :book: helped you :smiling_face_with_smiling_eyes: !"
                     )
    

    def changePath(self, newPath):
        """
        change the path to the complete path to a new working directory

        Parameters
        ----------
        newPath : str
            the new path.

        Returns
        -------
        str
            Print message indicating the new path.

        """
        self.path = newPath
        chdir(self.path)
        return print("path has changed to", self.path)

    @setDocStrAndApply  
    def setMatPltLibDefaults(self, default=False, style="seaborn-whitegrid",
                             rc_fig=dict(autolayout=True),
                             rc_ax=dict(labelweight="bold",
                                        labelsize="large",
                                        titleweight="bold",
                                        titlesize=14,
                                        titlepad=15,),
                             custom_cycler=(cycler(color=['c', 'm','y', 'k'])+
                                            cycler(lw=[1, 2, 3, 4])),
                             font = {'family' : 'sans-serif',
                                     'weight' : 'bold'}
                              ):

        matplotlib_params(default=default, style=style, rc_ax=rc_ax,
                          rc_fig=rc_fig, custom_cycler=custom_cycler,
                          font=font
                          )
        return print("plt default parameters have been setted")

    def plotBondAngleDist(self, filename='out.bad',
                          angLab=(r'S-Sb-S',  r'Sb-S-Sb'),
                          posAngLab=(0.01, .10, 0.10),
                          colors=('c', 'm'),
                          fig_param = dict(figsize=(6.4, 4.5),
                                           layout='constrained',
                                           edgecolor='white',
                                           facecolor='white',
                                           dpi=160, linewidth=0.0
                                           ),
                          plot_param = dict(marker=None, linestyle='-',
                                            linewidth=None, markersize=None),
                                  set_lim=True, limits=[(0,5), (0,)]
                        ):
        """
            Plot bond angles distribution.
    
        Parameters
        ----------
        filename: str, optional
            the name of the file.
            The default is 'out.bad'
        angLab : tuple, optional
            the different types of bonds. The default is (r'S-Sb-S',
                                                          r'Sb-S-Sb').
        posAngLab : tuple, optional
            tuple of size 3. First value is the horizontal position, 2nd is the
            vertical position and the 3rd the interline.
            The default is (0.01, 0.10, 0.10).
        colors : tuple, optional
            first value is the color of the fris bond plot and
            the second for the second. The default is ('c', 'm').
            Note :
                The color sequence must much cycler when the latter is setted
        fig_param : dict, optional
            Parameters controlling figure layout.
            The default is dict(figsize=(6.4, 4.5),
                                             layout='constrained',
                                             edgecolor='white',
                                             facecolor='white',
                                             dpi=160, linewidth=0.0
                                             ).
        plot_param : dict, optional
            parameters controlling the plot line.
            The default is dict(marker=None, linestyle='-',
                                             linewidth=None, markersize=None).
        set_lim : bool, optional
            If true, activate xlim, ylim. The default is True.
        limits : list, optional
            list of the tuple of xlim and ylim. The default is [(0,5), (0,)].
    
        Returns
        -------
        None.
    
        """
        angle, bad1, bad2, *_  = np.loadtxt(filename, unpack=True)
        
        # settin of figure and ax
        fig, ax = plt.subplots(**fig_param)
        
        # parameters of the labels of the axes
        ax_lab = (r"Angle (°)", r"")   
        
        # plotting of the data E_1 and E_2 as function of time
        # plotting of the axis
        ax.plot(angle, bad1, color=colors[0], **plot_param)
        ax.plot(angle, bad2, color=colors[1], **plot_param)
        
            
        
        # setting of the label to the axis
        ax.set_xlabel(ax_lab[0], fontsize=20)
        ax.set_ylabel(ax_lab[1], fontsize=20)
        
        # fill between plot
        ax.fill_between(angle, bad1,
                        where=bad1>=0.,alpha=.1, linewidth=0.1)
        ax.fill_between(angle, bad2,
                        where=bad2>=0.,alpha=.1, linewidth=0)
        
        #xticks
        ax.tick_params(labelcolor='black', labelsize=18, width=3)
        
        
        # Adding of additional informations on the plots
        i=0
        for label, color in zip(angLab, colors):
            ax.text(posAngLab[0], posAngLab[1] + i*posAngLab[2], angLab[0],
                    fontsize=18, color=color, transform=ax.transAxes)
            i +=1
        
        # saving figure and closing
        plt.autoscale(True,'both')
        ax.set_ylim(*limits[0])
        ax.set_xlim(*limits[1])
        #plt.yscale('log')
        #ax.set_xscale('logit')
        plt.savefig('bad.png')
        plt.show()
        plt.close(fig)

        
    def plotPartialgOfR(self, NbSpecies=2, filenames=('out.rdf',),
                        putFilLegend=False,
                        fileLab=(r'',), posText=(0.46, .40, 0.12), axe=1, 
                        bonds=('Sb-Sb', 'Sb-S', 'S-S'), posTextBonds=(.8, .8),
                        colors=('orange', 'C1', 'steelblue'), Textsize=18,
                        fig_param=dict(figsize=(5.4, 8), layout='constrained',
                        edgecolor='white', facecolor='white', dpi=200,
                        linewidth=3.5),
                        plot_param = dict(marker=None, linestyle=None,
                                          linewidth=2.0, markersize=None)
                        ):
        """
        

        Parameters
        ----------
        NbSpecies : int, optional
            The number of species. The default is 2.
        filenames : tuple of str, optional
            list of the filenames containing different g(r).
            The default is ['out.rdf'].
        putFilLegend : bool, optional
            Put on legending. The default is False.
        fileLab : tuple of str, optional
            The list of the labels to make the difference between the diffent 
            g(r). The default is [r''].
        posText : tuple, optional
            Respectively positions along x, and y, and the interline
            (associated to filelab)
            The default is (0.46, .40, 0.12).
        axe : int, optional
            the figure on which the labels will be write on.
            0 for the first (in row), 1 for the second ...
            The default is 1.
        bonds : tuple of str, optional
            the different types of bonds.
            The default is ('Sb-Sb', 'Sb-S', 'S-S').
        posTextBonds : tuple, optional
            Respectively positions along x, and y.
            (associated to bonds)
            The default is (0.01, .10).
        colors : tup of str, optional
            The color used for (resp.) the different type of bonds.
            The default is ('orange', 'C1', 'steelblue').
        Textsize : int, optional
            the size of these latters. The default is 18.
        fig_param : dict, optional
            control figure layout.
            The default is dict(figsize=(5.4, 8),
                                layout='constrained',
                                edgecolor='white', facecolor='white',
                                dpi=200, linewidth=3.5).
        plot_param : dict, optional
            control plot style.
            The default is dict(marker=None, linestyle=None,
                                linewidth=2.0, markersize=None).

        Returns
        -------
        None.

        """
        #############################################################
        # 22/12/2023 18h18 : color cycling to be optimized 
        # settings for the plot
        fig3d, ax = plt.subplots(nrows=len(bonds),
                                 ncols=1,
                                 **fig_param,
                                 sharex=True
                                 )

        #fig3d.subplots_adjust(bottom=0.15, left=0.2)
    
        #parameter for the labels of the different axis
        ax_lab = (r'r ($\AA$)', r'Partial pair correlation functions')
        n=0
        for file in filenames:
            data  = np.loadtxt(file, unpack=True)
            j=1
            for i in range(len(bonds)):
                ax[i].plot(data[0], data[j], label=bonds[i],
                           color=colors[n] ,**plot_param)
                j += 2
                ax[i].tick_params(labelcolor='black', labelsize=18, width=3)
                ax[i].text(*posTextBonds, bonds[i], fontsize=18,
                           transform=ax[i].transAxes)
            n+=1
                
        if putFilLegend:
            n=0
            for lab in fileLab:
                ax[axe].text(posText[0], posText[1] + n*posText[2],
                             lab, fontsize=Textsize, color=colors[n],
                             transform=ax[axe].transAxes)
                
                ax[axe].set_ylabel(ax_lab[1], fontsize=20)
                
                n+=1
            
        ax[-1].set_xlabel(ax_lab[0], fontsize=20)

        # other operations
        plt.autoscale(True,'both')
        plt.savefig('gdr.png')
        plt.show()
        plt.close(fig3d)

        
    def plotMSDAndDiffCoef(self, filenames=("msd.dat", "diff.dat"),
                           nspecies=2,
                           fig_param=dict(figsize=(6.4, 4.8),
                                             layout='constrained',
                                             edgecolor='white',
                                             facecolor='white',
                                             dpi=150, linewidth=2.0
                                             ),
                           tup_plot = (r"$S$",r"Sb"),
                           plot_param = dict(marker=None, linestyle=None,
                                             linewidth=3.0, markersize=None)
                       ):
        
        
        """
        plot the MSD and the diffusion coefficient as function of time.
        
        Parameters
        ----------
        filenames : tuple of str, optional
            list of size 2. first element is the name of the file containing
            MSD data. The second is for the diffusion coefficient.
            The default is ("msd.dat", "diff.dat").
        nspecies : int, optional
            the number of species composing the system. The default is 2.
        fig_param : dict, optional
            control figure layout.
            The default is dict(figsize=(6.4, 4.8),
                                layout='constrained',
                                edgecolor='white', facecolor='white',
                                dpi=150, linewidth=2.0 
                                ).
        tup_plot : tuple, optional
            The labels to differentiate the species (should respect the order
                                                     in the file).
            The default is (r"$S$",r"Sb").
        plot_param : dict, optional
            Control plot style. The default is dict(marker=None,
                                                    linestyle=None,
                                                    linewidth=3.0,
                                                    markersize=None).
    
        Returns
        -------
        None.
    
        """
        
        fig3d, ax = plt.subplots(nrows=2, ncols=1,  **fig_param,sharex=True)
        
        # parameters for the legend of the different axis
        dic_legend= dict(loc=0, shadow=False, fontsize=15)
    
        #parameter for the labels of the different axis
        ax_lab = ((u'', r'MSD ($\AA^2$)'),
                  (u'time (ps)', r'D ($cm^2.s^{-1}$)'))
        
        for i in range(2):
            data = np.loadtxt(filenames[i], unpack=True)
            
            for j in range(1, nspecies+1):
                ax[i].plot(data[0], data[j], label=tup_plot[j-1], **plot_param)
        
            # setting of the label to the axis
            ax[i].set_xlabel(ax_lab[i][0], fontsize=20)
            ax[i].set_ylabel(ax_lab[i][1], fontsize=20)
        
            #xticks
            ax[i].tick_params(labelcolor='black', labelsize=18, width=3)
        
            # setting of the legend
            ax[i].legend(**dic_legend)
    
        # other operations
        plt.autoscale(True,'both')
        plt.savefig('MSDandDiff.png')
        plt.show()
        plt.close(fig3d)

    
        
    def plotStress(self, dispMeanStress=True, lastConfig=200,
                   repTime=False, tstep=5.0, plotDiff=False, dpi=150
                   ):
        """
        plot the evolution of the stress

        Parameters
        ----------
        dispMeanStress : bool, optional
            If True the evolution of the average stress will also be plotted.
            The default is True.
        lastConfig : int, optional
            the meanPress calculated on the last lastConfig will be calculated.
            and printed.
            The default is 200.
        repTime : bool, optional
            If True, fix time accumulation activation problem.
            The default is false.
        tstep : float, optional 
            The time step in a.u..
        plotDiff : bool, optional
            If True abs(P - Pmean) will also be plotted on a new axe.
            The default is False.
        dpi : int, flaot, optional
            resolution of the image. The default is 150.
        Returns
        -------
        None.

        """
        tstep = tstep * 0.024 * 1e-3 # ps
        plotPressCelldim(plotMeanStress=dispMeanStress, lastConfig=200,
                         both=False, repTime=repTime, tstep=tstep,
                         plotDiff=plotDiff,dpi=dpi)
        plt.rcdefaults()

    def plotStressAndCellDim(self, plotMeanStress=True, repTime=False,
                             tstep = 5.0, dpi=150
                             ):
        """
        plot both the evolution of the stress and of the cell dimension.
        Useful in NPT simulation.

        Parameters
        ----------
        plotMeanStress : bool, optional
            If True the evolution of the average stress will also be plotted.
            The default is True.
        repTime : bool, optional
            If True, fix time accumulation activation problem.
            The default is false.
        tstep : float, optional 
            The time step in ps.
        dpi : int, flaot, optional
            resolution of the image. The default is 150.
        Returns
        -------
        None.

        """
        tstep = tstep * 0.024 * 1e-3 # ps
        plotPressCelldim(plotMeanStress=True, lastConfig=200, both=True,
                         repTime=repTime, tstep=tstep, dpi=dpi)

###############################################################################
#                                   EXCITON                                   #
###############################################################################
class Exciton:
    """
        Contains methods to plot exciton MSD
    """
    def __init__(self, path:str):
        """


        Parameters
        ----------
        path : str
            path to the directory of interest..

        Returns
        -------
        None.

        """
        self.path = path
        chdir(self.path)
        self.fileSpectr='transMoment.dat'
        print("[bold magenta][cyan]*[/][/]"*79,
              "[bold magenta][cyan]*[/][/]"+' '*35+"[bold magenta]EXCITON[/]"+
              ' '*35+"[bold magenta][cyan]*[/][/]",
              "[bold magenta][cyan]*[/][/]"*79,
              sep='\n')
        print(":Warning:")
        print(
            """

            For informations about all the available methods :
                -> Apply getInfos()
            For informations about the arguments of a given method :
                -> Apply infoMeth(obj.method); with obj the object you created
                    and its method.
                                        
                    Or simply use :
                        >> print(instance.method.__doc__)
            """
            )

    def infoMeth(self, method):
        """
        Print the documenation about a given method.

        Parameters
        ----------
        method : method
            The method of interest.
            e.g: if func() is a method of the class object h / h = Test()
            h.infoMeth(h.func) print the documentation of the method func()
        Returns
        -------
        None.

        """
        print(method.__doc__)

    def getInfos(self):
        msg ="""
        This Class Contains the following methods:
            changePath() :
                If you want to change the path from the one you
                defined when creating the object structuralProp().
            
            setMatPltLibDefaults(): 
                That set the layout of the plt figures.
                Note that after the end of each plot it is reset to the
                default matplotlib prop.
                
            MSDExciton() :
                plot the MSD of the Exciton.
                
            XtractransMom():
                Extract transition informations (LUMO, HOMO, TRANSITION MOMENT)
                from CPMD MATRIX.DAT file from different directories.
                
            lifetime():
                Calculate lifetime (if provided emission spectrum data) and/or 
                plot it.
            
            multiPlotSpectr():
                plot multiple transMoment.dat file (for emission and/or for 
                                                    absorption)
                
                
                
        Note :
            There is no method loadData since data are automatically
            loaded.
            
        """
        print(msg)
        return print("Hope this doc"+
                     " :book: helped you :smiling_face_with_smiling_eyes: !"
                     )

    def changePath(self, newPath):
        """
        change the path to the complete path to a new working directory

        Parameters
        ----------
        newPath : str
            the new path.

        Returns
        -------
        str
            Print message indicating the new path.

        """
        self.path = newPath
        chdir(self.path)
        return print(":Warning: path has changed to", self.path)

    @setDocStrAndApply  
    def setMatPltLibDefaults(self, default=False, style="seaborn-whitegrid",
                             rc_fig=dict(autolayout=True),
                             rc_ax=dict(labelweight="bold",
                                        labelsize="large",
                                        titleweight="bold",
                                        titlesize=14,
                                        titlepad=15,),
                             custom_cycler=(cycler(color=['c', 'm','y', 'k'])+
                                            cycler(lw=[1, 2, 3, 4])),
                             font = {'family' : 'sans-serif',
                                     'weight' : 'bold'}
                              ):

        matplotlib_params(default=default, style=style, rc_ax=rc_ax,
                          rc_fig=rc_fig, custom_cycler=custom_cycler,
                          font=font
                          )
        return print("plt default parameters have been setted")

    def MSDExciton(self, filenames=["msd_3d.time",
                              "msd_x.time",
                              "msd_y.time",
                              "msd_z.time"
                              ],
                   dimensions=(3, 1, 1, 1), ax_name=["3D", "X", "Y", "Z"],
                   nrows=2, ncols=2,
                   posText=[(.6, 0.01), (.6, .01),
                            (.6, .01), (.6, .01)
                            ],
                   fig_param = dict(figsize=(15., 8.0), layout='constrained',
                                    edgecolor='white', facecolor='white',
                                    dpi=160, linewidth=1.0
                                    ),
                   plot_param = dict(marker=None, linestyle=None,
                                     linewidth=3.0, markersize=None
                                     ),
                   sharex=True
                   ):
        """
        plot MSD of the Exciton

        Parameters
        ----------
        filenames : list or tuple of str, optional
            Liste of the name of the msd files to be plotted.
            The default is ["msd_3d.time",
                            "msd_x.time",
                            "msd_y.time",
                            "msd_z.time" 
                            ].
        dimensions : tuple of int, optional
            Dimensionality of the MSD with respect to the provided files.
            The default is (3, 1, 1, 1).
        ax_name : list or tuple of str, optional
            The names or labels of the different axes respectively to the given
            files (in uppercase). The default is ["3D", "X", "Y", "Z"].
        nrows : int, optional
            Number of axes in row. The default is 2.
        ncols : int, optional
            Number of axes in col. The default is 2.
        posText : list of tuple of folat, optional
            Position of the text writted on the figures.
            The default is [(.6, .01), (.6, .01),
                            (.6, .01), (.6, .01)
                          ].
        fig_param : dict, optional
            Figure features.
            The default is dict(figsize=(15., 8.0), layout='constrained',
                                edgecolor='white', facecolor='white', dpi=160,
                                linewidth=1.0
                                ).
        plot_param : dict, optional
            Plot features. The default is dict(marker=None, linestyle=None,
                                               linewidth=3.0, markersize=None
                                               ).
        sharex : bool, optional
            Indicate if x axes are shared. The default is True.

        Returns
        -------
        None.

        """
        print("\t[bold green underline]MEAN SQUARED DISPLACEMENT[/]\n")
    
        ps = 1e-12
        angs2Tocm2 = 1e-16 #cm**2
         
        # settings for the plot
        fig, ax = plt.subplots(nrows=nrows, ncols=ncols,
                               **fig_param,sharex=sharex
                               )
        _ax = []
        for axi in generatorAxeIndex(ax, ncols, nrows):
            _ax.append(axi)
        
        idx=0
        print("\t\tPlots :")
        
        # importing and defining quantities to plot the MSD
        for file, dim, axe in zip(filenames, dimensions, ax_name):
            
            time, msd  = np.loadtxt(file, unpack=True)
            print(f"\t\t\t[bold green]-> MSD {axe}[/]")

            # fittings
            popt, pcov = curve_fit(linear, time, msd)
            
            # mapping 
            fit = linear(time, popt[0])
    
            # calcul coefs
            D = popt[0] * angs2Tocm2 / (2. * dim * ps)
            
            
            # labels for the different plots of differents axis
            dic_plot={}
            dic_plot["3D"] = r"3D"
            dic_plot["X"] = r"along x"
            dic_plot["Y"] = r"along y"
            dic_plot["Z"] = r"along z"
            
            # parameters for the legend of the different axis
            dic_legend= dict(loc=0, shadow=False, fontsize=16)
        
            #parameter for the labels of the different axis
            ax_lab = (u'time (ps)', r'MSD ($\AA^2$)')
        
            _ax[idx].plot(time, msd, label=dic_plot[axe], **plot_param)
            _ax[idx].plot(time, fit, label=r'linear fit', **plot_param)
            
            #xticks
            _ax[idx].tick_params(labelcolor='black', labelsize=18, width=3)
            
            # setting of the legend
            _ax[idx].legend(**dic_legend)
        
            # Adding of additional informations on the plots
            _ax[idx].text(*posText[idx], r'D = %.2e $cm^2 s^{-1}$' %D,
                      fontsize=16, color='k', transform=_ax[idx].transAxes)
    
            # setting of the label to the axis
            _ax[idx].set_xlabel(ax_lab[0], fontsize=20)
            _ax[idx].set_ylabel(ax_lab[1], fontsize=20)
            
            idx += 1
        print("Job done ! :Thumbs_up:")
        # other operations
        plt.autoscale(True,'both')
        plt.savefig('MSD.png')
        plt.show()
        plt.close(fig)

        
    def XtractransMom(self, pathCalc=None, nameTargetDir='2',
                      nameTargetFile="MATRIX.DAT",
                      outFileName='transMoment.dat',
                      numCols=(5, 6, 10)
                      ):
        """
        Extract transition informations (LUMO, HOMO, TRANSITION MOMENT) from
        CPMD MATRIX.DAT file from different directories.

        Parameters
        ----------
        pathCalc : str, optional
            Complete path to the cpmd out file MATRIX.dat.
            If None the path current path is used.
            The default is './'.
        nameTargetDir : str, optional
            Name of the targeted directory. The default is '2'.
        nameTargetFile : str, optional
            The name of the target file. The default is "MATRIX.DAT".
        outFileName : str, optional
            Desired name for the output file. The default is 'transMoment.dat'.
        numCols : str, optional
            The columns to the grab from the file. The default is (5, 6, 10).

        Returns
        -------
        None.

        """

        self.fileSpectr=outFileName   
        if pathCalc is None:
            path=getcwd() #self.path

        else:
            path=pathCalc
        
        try:
            open(outFileName,'r')
        except OSError:
            parserTM(pathCalc=path, nameTargetDir=nameTargetDir,
                     nameTargetFile=nameTargetFile,
                     outFileName=outFileName, numCols=numCols)
        else:
            print("File [bold green]'transMoment.dat'[/] already exists.")
            print("Would you like to update it  ([bold green]yes[/]" +
            " [bold purple]or[/] [bold red]no[/]) :question:")
            decision = input(">> ")
            if decision.strip().upper() in ['YES', "OUI", "Y", '']:
                parserTM(pathCalc=path, nameTargetDir=nameTargetDir,
                         nameTargetFile=nameTargetFile,
                         outFileName=outFileName, numCols=numCols)
                print("\nThe file has been updated. :Thumbs_up:")
            else:
                print(":Warning: [bold yellow]"+
                      "You decided to keep existing file unchanged.")
        
    def lifetime(self, D=8.0e-3, T=300, err=0.,
                minE=1.6, maxE=2.0, plot=True, pltColor="limegreen",
                barColor="dodgerblue",oscStrength=True, aFile=None, 
                abarColor=None, apltColor=None, aminE=None, amaxE=None,
                npoints=1000, fontsize=7, barwidth=0.01, 
                figsize=(5., 2.5), dpi=800, file='transMoment.dat', locLeg=0,
                xlim=None, ylim=None, labels=('', ''), printLifeT=True,
                histNRJ=False
                ):
        """
        Calculate lifetime (if provided emission spectrum data) and/or plot it.
        Note
        ----
        All energy must be in eV
        
        Parameters
        ----------
        D : float, optional
            diffusion coefficient in cm**2/s. The default is 8.0e-3.
        T : float, optional
            Temperature in K. The default is 300.
        err : float, optional
               err = Egap_Exp -  E_00.
               
               Egap_Exp : experimental optical gap
               E_00 : position of the 0-0 transition pic (in emission spectrum)
                      before correction, of course.
        minE : float, optional
            minimum desired value for E (after correction)
            (important for the plot).
            The default is 1.6.
        maxE : float, optional
            max val (after correction). The default is 2.0.
        plot : bool, optional
            indicate whether we want to plot the spectrum or not.
            The default is True.
        pltColor : str, optional
            color of the line plot. The default is "limegreen".
        barColor : str, optional
            color of the barplot. The default is "dodgerblue".
        oscStrength : bool, optional
            If True the oscillator strength(s) will be plotted on the spectrum.
            The default is True.
        aFile : str, optional
            Complete path to the file to be superposed on already existing
            spectrum. It can be for example adding an absorption spectron upon
            an emission spectrum.
            The default is None.
        abarColor : str, optional
            Color of oscillator strength bar for the supplemental spectrum.
            The default is None.
        apltColor : str, optional
            Color for the new spectrum. The default is None.
        aminE : float, optional
            Its min energy value (after correction). The default is None.
        amaxE : float, optional
            Its max energy value (after correction). The default is None.
        npoints : float, optional
            number of desierd points for the gaussian broadening.
            The default is 1000.
        fontsize : int, optional
            control the size of the legend font. The default is 7.
        barwidth : float, optional
            control the width of the bars. The default is 0.01.
        figsize : tuple, optional
            size of the figure. The default is (5., 2.5).
        dpi : int, optional
            figure resolution. The default is 800.
        dpi : int, optional
            image resolution. The default is 800.
        file : str, optional
            The name of the targeted file (concern only the fisrt spectrum not
                                           the added).
            The default is 'transMoment.dat'.
        locLeg : int, optional
            location of the legend. The default is 0.
        xlim : tuple of float, optional
            x axis limits. The default is None.
        ylim : tuple of float, optional
            y axis limits. The default is None.
        labels : tuple of str, optional
            name to make difference between the two superposed spectra.
            The default is ('', '').

        Returns
        -------
        None.


        """
        msg=("The file [bold green]{}[/] [bold red]does not exists[/] "+
              "or [bold red]you have forgotten[/] to extract it."
            )

        print("\t[bold green underline]LIFETIME AND SPECTRUM[/]\n")

        try:
            # obf = file.split("\\")[-1]
            # print(file)
            open(file,'r')
        except OSError:
            print(":Warning:")
            print(msg.format(file))

            print("Please extract before and retry after "+
                  ":smiling_face_with_smiling_eyes:")
            print("Here is the method to do it: XtractransMom()")
        else:
            E_homo, E_lumo, tm = np.loadtxt(file, unpack=True)

            sup=not(aFile==None)
            if sup:
                aE_homo, aE_lumo, atm = np.loadtxt(aFile, unpack=True)
            else:
                aE_homo=aE_lumo=atm=None

            A21 = calcTau(E_homo, E_lumo, tm, D, T, err, minE, maxE, plot,
                          pltColor, barColor, npoints, fontsize, barwidth,
                          figsize, dpi, spltColor=apltColor,
                          sbarColor=abarColor, oscStrength=oscStrength,
                          superpose=sup, sminE=aminE, smaxE=amaxE,
                          sE_homo=aE_homo, sE_lumo=aE_lumo, stm=atm,
                          locLeg=locLeg, xlim=xlim, ylim=ylim, labels=labels,
                          printLifeT=printLifeT, histNRJ=histNRJ
                          )

            print("Job done :Thumbs_up: :Ok:")


    def multiPlotSpectr(self, eps, T=300, minE=1.6, maxE=2.0,
                        aminE=None, amaxE=None, npoints=1000,
                        fontsize=7, figsize=(5., 2.5), dpi=800,
                        eFile='transMoment.dat', ecolor=("b", "r"),
                        aFile=None, elw=2.0 , els='-',
                        alw=2.0, als='-', alpha_e=0.8, alpha_a=0.8,
                        locLeg=0, xlim=None, ylim=None, elabels='', alabels='',
                        ploTot_e=True, color_e='b', ploTot_a=True, 
                        annotate=False, expFile=None, expLab=None,
                        marker='o', markersize=2., expAlpha=2.,color='b',
                        vshift:str=0.7, fmt:str='png', finMaxVal:bool=False,
                        intervals:tuple=None, targ:int=None, toeV=True, shft=.0,
                        sep=',', dec='.'
                        ):
        
        """
        Plot multiple transMoment.dat on a single axis (for emission and/or for
                                                        absorption)

        Parameters
        ----------
        eps : float, optional
               eps = Egap_Exp -  E_00.
               
               Egap_Exp : experimental optical gap
               E_00 : position of the 0-0 transition pic (in emission spectrum)
                      before correction, of course.
        T : float, optional
            Temperature. Important for gaussian broadening. The default is 300.
        minE : float, optional
            min value for the absorption energie in eV
            (correction taken into account). The default is 1.6.
        maxE : float, optional
            min value for emission energie in eV 
            (correction taken into account). The default is 2.0.
        aminE : float, optional
            for the absoprtion spectra if there is. The default is None.
        amaxE : float, optional
            for the absoprtion spectra if there is. The default is None.
        npoints : int, optional
            number of points for the gaussian broadening. The default is 1000.
        fontsize : int, optional
            Size of the legend. The default is 7.
        figsize : tuple of float, optional
            figure size. The default is (5., 2.5).
        dpi : int, optional
            Image resolution. The default is 800.
        eFile : str or tup of str or list of str, optional
            names of the files for the emission spectra.
            
            Note:
                the path to a file can be precised in its name.
            The default is 'transMoment.dat'.
        aFile : str or tup of str or list of str, optional
            names of the files for the absorption spectra.
            
            Note:
                the path to a file can be precised in its name.
            The default is None.
        elw : float, tup of float or list of float, optional
            Line width of the plots for the emission spectra.
            The default is 2.0

        els : str or tup of str or list of str, optional
            Line style of the plots for the emission spectra.
            The default is '-'.
        alw : float, tup of float or list of float, optional
            Line width of the plots for the absorption spectra.
            The default is 2.0

        als : str or tup of str or list of str, optional
            Line style of the plots for the absorption spectra.
            The default is '-'.
        
        alpha_e : float, optional
            Transparency of the plots (for emission)
            The default is 0.8
        alpha_a : float, optional
            Transparency of the plots (for absorption)
            The default is 0.8
        locLeg : int or str, optional
            location of the legend. The default is 0.
        xlim : float, optional
            limit to the x-axis. The default is None.
        ylim : flaot, optional
            limit to the y-axis. The default is None.
        elabels : str or tup of str, or list of str, optional
            labels with respect to the eFile.
            Note :
                elabels must be of same type as eFile
            The default is ''.
        alabels : str or tup of str, or list of str, optional
            labels with respect to the eFile.
            Note :
                alabels must be of same type as aFile
            The default is ''.
        ploTot_e: bool, optional
            Indicates whether or not plot total spectrum (for emission)
            The default is True.
        color_e : TYPE, optional
            Color for the total spectra for emission.
            The default is 'b'.
        ploTot_a: bool, optional
            Indicates whether or not plot total spectrum (for absorption)
            The default is True.
        annotate : bool, optional
            Activate annotation. The default is False.
        expFile : str, list of str, or tuple of str, optional
            Experimental data files to plotted. The default is None.
        expLab : str, list of str, or tuple of str, optional
            The labels for the differents exp files. The default is None.
            
            Note: Must be of same type (and size if type is list or tuple) as
                  as expFile
        marker : str, list of str, or tuple of str, optional
            plt style for every exp data to be plotted. The default is 'o'.
        markersize : float, list of float, or tuple of float, optional
            the size of the exp plots. The default is 2..
        expAlpha : int or float, optional
            the transparency of the plot. The default is 2..
        vshift: float, optional
            vertical shift of the experimental data.
            The default is 0.7
        color : str, list of str, or tuple of str, optional
            colors for the different exp plots. The default is 'b'.
        toeV: bool, optional.
            If the exp spectr data to be converted into eV.
            The default is True.
        shft: float, optional.
            The shift to apply to the experimental spectra.
            The default is 0.
        Returns
        -------
            None
        """
        
        pltEmi = (eFile != None)
        pltAbs = (aFile != None)

        w = convKToeV * T


        def calcAndPlot(E_lumo, E_homo, tm, _minE, _maxE, lab, alpha, lw, ls,
                        color=None, finMax:bool=False, intrvl:tuple=None,
                        target:int=None
                        ):
            E, DeltaE, f12, _  = calcspctrDat(eps, E_lumo, E_homo, tm,
                                               histNRJ=False)
            # gaussian broadening of the spectra
            Enrg, g=gaussBroadn(Eigen=E, minVal=_minE, maxVal=_maxE,
                                npoints=npoints, height=f12, width=w)

            ax.plot(Enrg, g / max(g), alpha=alpha, linewidth=lw,
                    linestyle=ls, label=lab, color=color)

            if finMax:
                maxE, maxg = finMaxNrj(array=np.array((Enrg,g)),
                                       inf=intrvl[0], sup=intrvl[1],
                                       target=target, col=0
                                       )
                ax.scatter(maxE, maxg/ max(g), color='k', alpha=0.7)#,linestyle=':')
                print("Max val in  [",*intrvl,f"] is {maxE}")
        
        if pltAbs or pltEmi:
            fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figsize, dpi=dpi)
        
            if pltEmi:
                homo, lumo, trMom = np.array([[],[], []])

                if isinstance(eFile, str):
                    eFile = eFile,
                    elabels = elabels,
                    ecolor = ecolor,
                    #intervals = intervals,
                    
                for file, lab, col in zip(eFile, elabels, ecolor):#,
                                            #intrvl):
                    E_homo, E_lumo, tm = np.loadtxt(file, unpack=True)
                    
                    homo = np.concatenate((homo, E_homo))
                    lumo = np.concatenate((lumo, E_lumo))
                    trMom = np.concatenate((trMom, tm))

                    calcAndPlot(E_lumo, E_homo, tm, minE, maxE, lab, alpha_e,
                                elw, els, color=col, finMax=finMaxVal,
                                intrvl=intervals, target=targ)

                if ploTot_e:     
                    calcAndPlot(lumo, homo, trMom, minE, maxE, "Total",
                                alpha=1.0, lw=2.0, ls='-', color=color_e)

            if pltAbs:
                homo, lumo, trMom = np.array([[],[], []])

                if isinstance(aFile, str):
                    aFile = aFile,
                    alabels = alabels,
                    ecolor = ecolor,

                for file, lab, col in zip(aFile, alabels, ecolor):
                    E_homo, E_lumo, tm = np.loadtxt(file, unpack=True)

                    if ploTot_a:
                        homo = np.concatenate((homo, E_homo))
                        lumo = np.concatenate((lumo, E_lumo))
                        trMom = np.concatenate((trMom, tm))

                    calcAndPlot(E_lumo, E_homo, tm, aminE, amaxE, lab, alpha_a,
                                alw, als, color=col)

                if ploTot_a:
                    calcAndPlot(lumo, homo, trMom, aminE, amaxE, "Total",
                                alpha=1., lw=2.0, ls='-')


            ax.set_xlabel('$E$ (eV)')
            ax.set_ylabel(r'Intensity (a. u.)')
            ax.tick_params(axis='both', direction='in',
                           length=4.0, width=1.0)

            print("The intensities are normalized regarding the "+
                  "highest peak in the spectrum.")

            if xlim==None and ylim==None:
                plt.autoscale(True,'both')
            elif not xlim==None:
                ax.set_xlim(xlim)
            elif not ylim==None:
                ax.set_ylim(ylim)
            else:
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)

            if annotate:
                # Annotation
                ax.annotate('0-0', xy=(1.88, 1.5), xytext=(1.855, 1.7),
                             arrowprops=dict(facecolor='black',
                                             arrowstyle='simple',
                                             shrinkA=3,
                                             shrinkB=5,
                                             mutation_scale=2.5)
                             )
                ax.annotate('0-1', xy=(1.69, 1.75), xytext=(1.665, 1.95),
                             arrowprops=dict(facecolor='black',
                                             arrowstyle='simple',
                                             shrinkA=3,
                                             shrinkB=5,
                                             mutation_scale=2.5)
                             )
                ax.annotate('0-2', xy=(1.52, 1.4), xytext=(1.495, 1.6),
                             arrowprops=dict(facecolor='black',
                                             arrowstyle='simple',
                                             shrinkA=3,
                                             shrinkB=5,
                                             mutation_scale=2.5))
            if expFile != None:
                if isinstance(expFile, (list, tuple)):
                    genFnameLab = zip(expFile, expLab, marker,
                                      markersize, expAlpha, color
                               )
                elif isinstance(expFile, str):
                    expFile = expFile,
                    expLab = expLab,
                    marker = marker,
                    markersize = markersize,
                    expAlpha = expAlpha,
                    color=color,
                    
                    genFnameLab = zip(expFile, expLab, marker,
                                      markersize, expAlpha, color)
                    
                for fname, lab, mk, mks, pha, col in genFnameLab:
                    data = pd.read_csv(fname, sep=sep, decimal=dec)
                    if toeV:
                        nrj = (h * c) / (data.iloc[:,0] * eV * 1e-9)
                        mag = data.iloc[:,1]/data.iloc[:,1].max()
                    else:
                        nrj = data.iloc[:,0] 
                        mag = data.iloc[:,1]/data.iloc[:,1].max()
                    
                    ax.plot(nrj, shft + mag, alpha=pha, linewidth=mks,
                            linestyle=mk, label=lab, color=col
                            )
                    if finMaxVal:
                        print(f"file {lab}")
                        maxE, maxg = finMaxNrj(array=np.array((nrj,shft + mag)),
                                               inf=intervals[0],
                                               sup=intervals[1],
                                               target=targ, col=0
                                               )
                        ax.scatter(maxE, maxg, color='k', alpha=0.7)
                        print("Max val in  [",*intervals,f"] is {maxE}")
                    # ax.scatter(dat.x, .7 + dat.y / max(dat.y), alpha=pha,
                    #            label=lab, marker=mk, s=mks
                    #            )
                
                
            #plt.xticks(np.arange(1.2, 2.1, 0.1))
            #plt.yticks([])
            fig.subplots_adjust(bottom=0.2, left=0.5)
            # for pos in [ "right", "top"]:
            #     ax.spines[pos].set_color('white')
                
            leg=dict(loc=locLeg, shadow=False, fontsize=fontsize)    
            plt.legend(**leg)
            
            plt.savefig(f"multISpctr.{fmt}", dpi=dpi, format=fmt)
            plt.show()
            plt.close(fig)

class FixGeo():
    """
        Class that try to reconstruct folded geometry.
        That is geometry where are structure converved atoms and some 
        isolated clusters due to the PBC
    """
    def __init__(self, distEps:dict=None, geofile:str="GEOMETRY.xyz",  
                 inSet:set = {7, 5, 16, 18,54, 41, 43, 52},
                 vec1:list=[26.0140, 0.0000, 0.0000],
                 vec2:list=[-6.7765, 61.4499, 0.0000],
                 vec3:list=[0.0000, 0.0000, 29.8766],
                 traj:bool=False,
                 path:str=None
                 ):
        
        """
        

        Parameters
        ----------
        distEps : dict
            The different species and corresponding characteristic interatomic
            distance disreagarding the specy of the with which the bonds are
            made.
            The default is    
            distEps = {}
            distEps['S'] = 1.95
            distEps['O'] = 1.8
            distEps['N'] = 1.8
            distEps['C'] = 1.9
            distEps['H'] = 1.6
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
        if distEps is None:
            distEps = {}
            distEps['S'] = 1.95
            distEps['O'] = 1.8
            distEps['N'] = 1.8
            distEps['C'] = 1.9
            distEps['H'] = 1.6

        self.distEps=distEps

        # Read geo file here when the input file is not for a trajectory
        # and initialise some characteristic of the system
        if path != None:
            chdir(path)
            
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
            geo = self.geoOut
        else:
            geo = self.geo
        dist_ij = distance(geo - geo[:, ind, np.newaxis])
        coordAtm = np.where(dist_ij<=self.distEps[self.specy[ind]])
        coordAtm = np.setdiff1d(coordAtm, ind)

        self.coordAtm[ind] += list(coordAtm)
        indInclus = (len(set(coordAtm).intersection(inclus)) != 0)
              
        if indInclus:
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
            np.savetxt('inSet.txt', list(self.inSet))
        
        self.saveGeo(geo, self.inSet, fileOut='GEOMETRY_inSet.xyz')

    def coordAndIndInc(self, ind, dist_ijOut):
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

    def addInclSuppExc(self, ind, geoOut, bigGeo):
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

    def fixTraj(self, fileTraj:str="TRAJEC.xyz"):
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
        """
        find indexes of atoms present in both a geometry of reference in
        in wich sidechains have been removed (for example), and an the 
        actual geometry

        Parameters
        ----------
        geoRef : np.array
            geometry of reference.
        eps : float, optional
            the typical interatomic distance. The default is 1.

        Returns
        -------
        indices : list
            list of indexes of atoms.

        """
        
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
                coordAtm = np.where(dist_ij<=self.distEps['C'])
                vois=list(coordAtm[0])
                if len(vois) == 2:
                    #print(coordAtm[0])
                    self.tobeSat.append([vois, i, round(dist_ij[vois[1]],2),
                                         self.specyUn[vois[1]]]
                                        )
  
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
            
        lenGeoUn = len(self.indices)# self.geoUn.shape[1]
        
  
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
  
                #saturate unsaturated C atoms 
                self.satur()
                
                for i in range(lenGeoUn):
                    file.write(data.format(self.specyUn[i], self.geoUn[0, i],
                                           self.geoUn[1, i], self.geoUn[2, i])
                               )
                for i in range(lenSatH):
                    a = self.satH[i]
                    file.write(data.format("H", a[0], a[1], a[2]))
        file.close()

        #read existing index
        # if self.indices is None:
        #     self.indices = np.loadtxt("indices.txt").astype('int64')
        # if self.satH is None:
        #     lenSatH = np.loadtxt("lenSatH").astype('int64')
        # else:
        #     lenSatH = len(self.satH)
        # if self.specyUn is None:  
        #     self.specyUn =[]
        #     for i in range(lenGeoUn):
        #         self.specyUn.append(self.specy[self.indices[i]]) 