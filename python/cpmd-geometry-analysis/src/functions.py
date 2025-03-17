# -*- coding: utf-8 -*-

"""
Created on Thu Nov  2 16:38:31 2023

@author: codiarra
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os 
from shutil import copytree
from cycler import cycler
from functools import wraps
from rich import print
import scipy.constants as cons
#from scipy.special import expit
from pandas import DataFrame

#______________________________________________________________________________
#                           IMPORTANT CONSTANTS
kB = cons.Boltzmann # J/K
hartree = 2.0 * 13.6 # eV
bohr = 5.2912*1e-11 # m
q = cons.e # C
hbar = cons.hbar # J.s
eps0 = cons.epsilon_0 #  (c/s)**2 s4 kg−1 m−3
c = cons.c # m/s
me =cons.m_e # kg
convKToeV =  kB / q # kB / eV (eV/K) conv from K to eV
#______________________________________________________________________________
#                       IMPORTANT GENERAL FUNCTIONS 
def repair_time(time, forStrs=False):
    """
    Fix time problem when forget to activate accumulators in CPMD simulations.

    Parameters
    ----------
    time : np.array
        Input raw time.
    forStrs : bool, optional
        True if for plot of Stress of celldim (in NPT for example).
        The default is False.

    Returns
    -------
    time
        monotonic time.

    """

    time_cor = 0
    time_out = []
    time_out.append(time[0])
    
    null = lambda *arg : 0 
    corrStr={True:corrTime, False:null}
    
    for t in range(1, len(time)):
        a = time[t-1]
        b=time[t]
        if (b <= a):
            time_cor = time_cor + a + corrStr[forStrs](b, time[t+1])

        time_out.append(b + time_cor)   
    return np.array(time_out[:])

def matplotlib_params(default=False, style="seaborn-whitegrid",
                      rc_fig=dict(autolayout=True),
                      rc_ax=dict(labelweight="bold",
                                 labelsize="large",
                                 titleweight="bold",
                                 titlesize=14,
                                 titlepad=15,),
                      custom_cycler=(cycler(color=['c', 'm','y', 'k']) +
                                     cycler(lw=[1, 2, 3, 4])),
                      font = {'family' : 'sans-serif',
                              'weight' : 'bold'
                              }
                      ):
        """
        set Matplotlib defaults

        Parameters
        ----------
        default : bool, optional
            If true the defaults parameter are directly used.
            Else, the parameter specified in argument will be take into 
            account.
            The default is False.
        style : str, optional
            Define the style of the figure layout (see the website).
            To get the list of available style, execute the following: 
                plt.style.available
            The default is "seaborn-whitegrid".
        rc_fig : dict, optional
            parameters of the figure. The default is dict(autolayout=True).
        rc_ax : dict, optional
            parameters of the axes. The default is dict(labelweight="bold",
                                                        labelsize="large",
                                                        titleweight="bold",
                                                        titlesize=14,
                                                        titlepad=10,).
        custom_cycler : cycler.cylcer, optional
            Can also take the value 'default', which give the default
            matplotlib color cycle.
            color cycling. The default is (cycler(color=['c', 'm','y', 'k']) + 
                                           cycler(lw=[1, 2, 3, 4])).

        Returns
        -------
        None.

        """
        plt.rcdefaults()

        if default:
            plt.rcdefaults()
        else:
            # Set Matplotlib defaults
            plt.style.use(style)
            plt.rc("figure", **rc_fig)
            if custom_cycler == 'default':
                custom_cycler = (cycler(color=['r', 'g', 'b', 'y']) +
                                  cycler(lw=[1, 2, 3, 4]))
            plt.rc("axes",prop_cycle=custom_cycler, **rc_ax)
            plt.rc('lines', linewidth=2, color='r')
            plt.rc('font', **font)  # pass in the font dict as kwargs

def setDocStrAndApply(function):
    docs = """
    set Matplotlib defaults

    Parameters
    ----------
    default : bool, optional
        If true the defaults parameter are directly used.
        Else, the parameter specified in argument will be take into 
        account.
        The default is False.
    style : str, optional
        Define the style of the figure layout (see the website).
        To get the list of available style, execute the following: 
            plt.style.available
        The default is "seaborn-whitegrid".
    rc_fig : dict, optional
        parameters of the figure. The default is dict(autolayout=True).
    rc_ax : dict, optional
        parameters of the axes. The default is dict(labelweight="bold",
                                                    labelsize="large",
                                                    titleweight="bold",
                                                    titlesize=14,
                                                    titlepad=10,).
    custom_cycler : cycler.cylcer, optional
        Can also take the value 'default', which give the default
        matplotlib color cycle.
        color cycling. The default is (cycler(color=['c', 'm','y', 'k']) + 
                                        cycler(lw=[1, 2, 3, 4])).

    Returns
    -------
    None.

    """
    function.__doc__ = docs
    return function


def listOfStrToListOfFloat(data, liste):
    
    """
    Convert list of element of type string convertible to float type into
    list of float.

    Parameters
    ----------
    data : list
        list of string elements.
    liste : list
        the output liste.

    Returns
    -------
    None.

    """
    
    for i in range(len(liste)):
        data[i].append(float(liste[i]))
        
def corrTime(time1, time2):
        return time2 - 2*time1

def delColSkip(line, skip):
    """
    delete given columns from the file

    Parameters
    ----------
    line : list
        list of one row of the data.
    skip : int, list or tuple
        columns to be deleted.

    Returns
    -------
    None.

    """       
        
    itr=0 # iterant
    if type(skip) is list or type(skip) is tuple:
        for number in skip:
            del line[number - itr]
            itr += 1
    elif type(skip) is int:
            del line[skip]
    
def grabCol(line, grepColumn, coldtype):
        return coldtype(line[grepColumn])
  
def importData(filename, ncols, com="#", ignore=None, sep=' ', skip=None,
               repTime=False, grepColumn=None, coldtype=str):
    """
    Import data  in a variable of type array.
    An home-made alternative to np.loadtxt

    Parameters
    ----------
    filename : string
        Name of the file.
    ncols : int
        The total number of columns in data file.
    com : string, optional
        The comments symbol. Indicate line to be ignired.
        The default is "#".
    ignore : int, tuple or list, or str, optional
        In case int, tuple, or list: The numbers of the lines to be ignored.
        In case str the format is 'debut:fin'.
        The default is None.
    sep : string, optional
        Columns delimiter character.
        The default is ' '.
    skip : int, list or tuple, optional
        If specified the column number(s)  skip will be ignored.
        for e.g.: skip=0 then the first column will be ignored.
                  skip=[0, 1] or (0, 1) the two first columns we be ignored.           
        The default is None.
    repTime : bool, optional
        If time (accumulator problem during MD simulation) to be repaired.
        The default is False.
    grepColumn : int, optional
        Will keep if specified the column number indicated keepAppart
        in list that will be returned. The default is None.
    coldtype : type, optional
        type of the column. Can only be either str, float or int.
        The default is str.
    

    Returns
    -------
    np.array
        array of the data imported.

    """
    n = 0 # iterant for the lines to be ignored

    # readjust the the numbers of cols cause when some will be dropped
    if type(skip) is int:
        ncols = ncols - 1
    elif type(skip) in (tuple, list):
        ncols = ncols - len(skip)

    # creation and initialisation of the list data
    data_in = []
    for i in range(ncols):
        data_in.append([])

    try :
        open(filename, 'r')
    except :
        print(f":Warning: {filename} does not exist !")
    else:
        of = open(filename, 'r')
        print(f"\t[bold green]READING {filename}[/]")
    
    iter = -1 # iterant for the number of the last line red

    grep = [] # only important when one to grep one cols dropped and return it

    # restructure ignore data
    if isinstance(ignore, int):
        ignore = ignore,
    elif type(ignore) is str:
        splIgnore = ignore.split(':')
        debut, fin = splIgnore
        ignore = []
        for i in range(int(debut), int(fin)+1):
            ignore = ignore.append(i)

    # operations really begin from here
    while True:
        line = of.readline()
        iter += 1
        if line:
            if (com in line):
                # commented lines are ignored
                line = of.readline()
                iter += 1
            elif (iter == ignore[n]):
                # some specified lines are also ignored
                line = of.readline()
                iter += 1
                if n+1 < len(ignore):
                    n += 1
            else :
                line = line.replace('D','E')  # e.g.:  1D08 to 1E08
                line = line.replace('\n','')
                line = line.split(sep)
                
                te = line[0] # cause don't found this special symbol to be rm
                while te in line:
                    line.remove(te)

                if type(grepColumn) != type(None):
                    # grep specified col
                    grep.append(grabCol(line, grepColumn, coldtype))

                if type(skip) != type(None):
                    # del specified cols
                    delColSkip(line, skip)

                # Structure data, i.e.: str -> float
                listOfStrToListOfFloat(data_in, line) # conv line to float type
        else:
            break
        
    # fix time accumulation problem
    if repTime:
        data_in[0] = repair_time(data_in[0])  
     
    data_in = np.array(data_in)
    
    if grepColumn is not None:
        return grep, data_in[:]
    else: 
        return data_in[:]

                 
#______________________________________________________________________________
#                              energies plot                                  #
#______________________________________________________________________________

def _plots(data, ax, ax_lab, legend_dic, lab_param, plot_param):
    
    """
    plot of the energies

    Parameters
    ----------
    data : tuple
        argument, i.e, data to plot on the same axes.
    ax : matplotlib.axes._axes.Axes
        DESCRIPTION.
    ax_lab : tup
        tuple of the label of ax.
    legend_dic : dict
        parameter of the legending.
    lab_param : dic
        parameters for the plot.
    plot_param : dict
        parameters of the plots.

    Returns
    -------
    None.

    """
        
    for i in range(1, len(data)):
        # plotting of the axis
        ax.plot(data[0], data[i], label=lab_param[i-1], **plot_param)
        # setting of the label to the axis
        ax.set_xlabel(ax_lab[0], fontsize=20)
        ax.set_ylabel(ax_lab[1], fontsize=20)
        #xticks
        ax.tick_params(labelcolor='black', labelsize=18, width=3)
        # setting of the legend
        ax.legend(**legend_dic)
        
conv_kin = lambda N_ions : (3./2.) * (N_ions) * (kB/q) / hartree # K to hartree
conv_fic = lambda N_e : (2.0 * (hartree * q)) / (N_e * kB) # E_fic to K 
        
def _plot2(data, ax, ax_lab, legend_dic, lab_param, plot_param, keys):
    
    """
    plot of the energies

    Parameters
    ----------
    data : tuple
        argument, i.e, data to plot on the same axes.
    ax : matplotlib.axes._axes.Axes
        DESCRIPTION.
    ax_lab : tup
        tuple of the label of ax.
    legend_dic : dict
        parameter of the legending.
    lab_param : dic
        parameters for the plot.
    plot_param : dict
        parameters of the plots.

    Returns
    -------
    None.

    """
        
    for i in keys:
        # plotting of the axis
        ax.plot(data[0], data[i], label=lab_param[i], **plot_param)
        # setting of the label to the axis
        ax.set_xlabel(ax_lab[0], fontsize=20)
        ax.set_ylabel(ax_lab[1], fontsize=20)
        #xticks
        ax.tick_params(labelcolor='black', labelsize=18, width=3)
        # setting of the legend
        ax.legend(**legend_dic)

def structInpData2(data, NPT=False):
    """
    Structures the data for NPT or NVE plot

    Parameters
    ----------
    data : np.array
        input data (directly the data in ENERGIES.
    NPT : bool, optional
        Indicates wether or not its for a NPT calculation. The default is True.

    Returns
    -------
    data_out : pandas.DataFrame
        Strcutured data.

    """
    dic_colsnames = {}
    dic_colsnames[False] = ("time", "EkinFic", "Tions",
                            "Epot", "Eclas", "Econs","_")
    dic_colsnames[True] = ("time", "EkinFic", "EkinCell", "Tions",
                            "Epot", "Eclas", "Econs", '_')
    data_out = DataFrame()
    
    for colIndx in range(data.shape[0]-1):
        data_out[dic_colsnames[NPT][colIndx]]=data[colIndx]
    
    return data_out

def setPlotLabels(NPT=False):
    
    dic_plot={} # contains labels of corresponding plot 
    plot_comp={}
    
    if NPT:
          dic_plot['EkinCell'] =r'E$_{Kin, Cell}$'
          
          plot_comp['Fict'] = ("EkinIons", "EkinFic", "EkinCell")
    else:
        dic_plot["Fict"] = (r"$E_{Kin, ions}$", r"$E_{Kin, fic}$")
        plot_comp['Fict'] = ("EkinIons", "EkinFic")
    
    # labels for the different plots of differents axis
    dic_plot['EkinIons'] = r"$E_{Kin, ions}$"
    dic_plot['EkinFic'] =r"$E_{Kin, fic}$"
    
    dic_plot["Econs"] = r"$E_{cons}$"
    dic_plot["Eclas"] = r"$E_{clas}$"
    dic_plot["Epot"] = r"$E_{pot}$"
    dic_plot["DeltaEcons"] = r"$\frac{\Delta E_{cons}}{E_{cons}(0)}$"
    dic_plot["Tions"] = r"$T_{ions}$"
    dic_plot["Tfic"] = r"$T_{fic}$"
    
    plot_comp['Cons'] = ("Econs", "Eclas", "Epot")
    plot_comp['Vari'] = ("DeltaEcons",)
    plot_comp['Temp'] = ("Tions", "Tfic")
    

    return dic_plot, plot_comp
        

def structInpData(data, N_ions, N_e, tstep=1., NPT=False):
    """
    Structures the data for NPT or NVE plot

    Parameters
    ----------
    data : np.array
        input data (directly the data in ENERGIES.
    N_ions : int
        Total number of ions.
    N_e : int
        total number of electrons.
    NPT : bool, optional
        Indicates wether or not its for a NPT calculation. The default is True.
    tstep : int, optional
        The time step in a.u..
        The default is 1.

    Returns
    -------
    data_dic : dict
        dictionary of the data structure and ready to be used for plot.
    dic_plot : dict
        dictionary of the labels of plots.
    meanTions : np.array
        mean temperature (from t=0 to the current time).

    """
    
    data_dic = {} # contains data restructured
    dic_plot={} # contains labels of corresponding plot 
    
    # Parameter for conversion
    convkin = (3./2.) * (N_ions) * (kB/q) / hartree # K to hartree
    convfic = (2.0 * (hartree * q)) / (N_e * kB) # E_fic to K 
    
    i=0 # if NPT i=1 else i=0 
        # cause in NPT there is intermediate cols between E_fic and T_ions
    
    T_fic = data[1] * convfic # fictitious temperature
    data[0] = data[0] * tstep  # time in ps 

    if NPT:
          i=1
          E_kin = convkin * data[3] # E_kin,ions in hartree
          meanTions = _meanT(data[3]) # compting mean temperature
          data_dic[1,0] = (data[0], E_kin, data[1], data[2])
          dic_plot[1,0] = (r"$E_{Kin, ions}$",
                          r"$E_{Kin, fic}$",
                          r'E$_{Kin, Cell}$')
    else:
        E_kin = conv_kin * data[2]
        meanTions = _meanT(data[2])
        data_dic[1,0] = (data[0], E_kin, data[1])
        dic_plot[1,0] = (r"$E_{Kin, ions}$", r"$E_{Kin, fic}$")
         
    data_dic[0,0] = (data[0], data[5+i], data[4+i], data[3+i]) #multiplot
    delta_Econs = (data[5+i] - data[5+i][0]) / data[5+i][0]
    data_dic[0,1] = (data[0], delta_Econs) #simple plot
    data_dic[1,1] = (data[0], data[2+i], T_fic)
    
    # labels for the different plots of differents axis
    dic_plot[0,0] = (r"$E_{cons}$", r"$E_{clas}$", r"$E_{pot}$")
    dic_plot[0,1] = (r"$\frac{\Delta E_{cons}}{E_{cons}(0)}$",)
    dic_plot[1,1] = (r"$T_{ions}$",r"$T_{fic}$")
    
    
    
    
    return data_dic, dic_plot, meanTions



#______________________________________________________________________________
#                       Exploitation des datas NPT                            #  
#  __________________________________________________________________________ #
#   Originally writted By Achille Lambrecht (some modification was brought to #
#    it by myself                                                             #  
#                                                                             #  
#______________________________________________________________________________

def readSTRESS(time, name):
    """
    Read the Stress components in the STRESS file from CPMD
    ------
    Param 
    ------
    time : list
        the time (will be modified in output)
    name :str
        path to file
    ------
    Return
    ------ 
        list (Stress Components)
    """
    with open(name,'r',encoding='utf8') as f:
        i=0
        STRESS=[]
        for line in f:
            line=line.split()
            if line[0]=='TOTAL':
                time.append(int(line[-1])) 
            else:
                if i<=2:
                    s=float(line[i])
                    STRESS.append(s)
                    i+=1
                else :
                    i=0
                    s=float(line[i])
                    STRESS.append(s)
                    i+=1
    return STRESS                
    
def CompPressure(allcomp):
    """
    Compute the pressure from the stress tensor components
    ------
    Param
    ------
    list (Stress Components)
    ------
    Return
    ------
    list (Pressure)
    """
    j = 1
    p = 0.0
    pressure = []
    for i in range(len(allcomp)):
        p += allcomp[i]
        j += 1
        if j == 3:
            p = (p+allcomp[i])/3
            pressure.append(p)
            j = 0
            p = 0.0
    return pressure

def Dimcell(name):
    """
    Extract the size of the a direction of the simulation cell
    -----
    Param
    -----
    str (path to file)
    -----
    Return
    -----
    list (dim cell)
    """
    with open(name, 'r', encoding='utf8') as f:
        alat = []
        time = []
        for line in f:
            line = line.split()
            if (line[0] == 'CELL'):
                time.append(int(line[-1]))
            elif (line[0] == '0.000000'):
                continue
            else:
                a=float(line[0])
                alat.append(a)
    return time, alat
            

def plotBoth(*arg, plotMeanStress=True, lastConfig=200, repTime=False, tstep,
             dpi
             ):
    """
    function that plot both the stress and the cell dimension as function of 
    the time. Will print the averaged stress on the last  lastConfig

    Parameters
    ----------
    *arg : arg
        arguments.
    plotMeanStress : bool, optional
        If True the evolution of the averaged stress will be plotted.
        The default is True.
    lastConfig : int, optional
        number of last pressure time point on which we want to average.
        The default is 200.
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
    
    timecell, dim = Dimcell(r'CELL')
     #timecell=np.linspace(10*5*0.024/1000, len(dim)*5*10*0.024/1000,
     #                    len(dim))
    if repTime:
        timecell = repair_time(time=timecell, forStrs=True)
    
    timecell = np.array(timecell) * tstep
        
    # figure parameters
    fig_param = dict(figsize=(6.4, 5), layout='constrained',
                     edgecolor='white', facecolor='white',
                     dpi=dpi, linewidth=3.5)

    # settings for the plot
    fig3d, ax = plt.subplots(nrows=2, ncols=1,**fig_param, sharex=True)
    #fig3d.subplots_adjust(bottom=0.15, left=0.2)

    # parameter for plots
    plot_param = dict(marker=None, linestyle=None,
                      linewidth=1.7, markersize=None)
    

    #parameter for the labels of the different axis
    ax_lab = (r'Time (ps)', r'Pressure (kbar)', r'alat (A.U.)')
    

    ax[0].plot(arg[0], arg[1], label='Pressure',
               **plot_param, color='midnightblue')
    if plotMeanStress:
        meanPress = _meanP(arg[1])
        ax[0].plot(arg[0], meanPress, label="mean pressure",
                   **plot_param, color='cyan')
    ax[0].legend()
    
    ax[1].plot(timecell, dim, **plot_param, color='midnightblue')

    # setting of the label to the axis
    ax[0].set_ylabel(ax_lab[1], fontsize=20)
    ax[1].set_xlabel(ax_lab[0], fontsize=20)
    ax[1].set_ylabel(ax_lab[2], fontsize=20)

    #xticks
    ax[0].tick_params(labelcolor='black', labelsize=18, width=3)
    ax[1].tick_params(labelcolor='black', labelsize=18, width=3)

    plt.autoscale(True,'both')
    plt.savefig('NPTplot.png')
    plt.show()
    plt.close(fig3d)
    
    
def onlyPress(time, Pressure, plotMeanStress=True, lastConfig=200,
              plotDiff=False, dpi=150):
    """
    Function that plot only the stress as function of the time.
    Will print the averaged stress on the last  lastConfig

    Parameters
    ----------
    time : np.array
        time
    Pressure : np.array
        Pressure data.
    plotMeanStress : bool, optional
        If True the evolution of the averaged stress will be plotted.
        The default is True.
    lastConfig : int, optional
        number of last pressure time point on which we want to average.
        The default is 200.
    plotDiff : bool, optional
        If True abs(P - Pmean) will also be plotted on a new axe.
        The default is False.
    dpi : int, optional,
        fix the image resolution.
        The default is 150.
    Returns
    -------
    None.

    """

    # settings for the plot
    if plotDiff:
        # figure parameters
        fig_param = dict(figsize=(6.4, 5), layout='constrained',
                         edgecolor='white', facecolor='white',
                         dpi=dpi, linewidth=3.5)
        fig3d, ax = plt.subplots(nrows=2, ncols=1,**fig_param, sharex=True)
        _ax = ax[0]
    else:
        # figure parameters
        fig_param = dict(figsize=(6.4, 2.5), layout='constrained',
                         edgecolor='white', facecolor='white',
                         dpi=150, linewidth=3.5)
        fig3d, ax = plt.subplots(nrows=1, ncols=1,**fig_param)
        _ax = ax

    # parameter for plots
    plot_param = dict(marker=None, linestyle=None,
                      linewidth=1.7, markersize=None)
    
    #parameter for the labels of the different axis
    ax_lab = (r'Time (ps)', r'Pressure (kbar)', r'|P - Pmean|')
    

    _ax.plot(time, Pressure, label='Pressure',
             **plot_param, color='midnightblue')
    
    meanPress = _meanP(Pressure)

    if plotMeanStress:
        _ax.plot(time, meanPress, label="mean pressure",
                   **plot_param, color='cyan')
    _ax.legend()
    
    # setting of the label to the axis
    _ax.set_ylabel(ax_lab[1], fontsize=20)


    #xticks
    _ax.tick_params(labelcolor='black', labelsize=16, width=3)
    if plotDiff:
        ax[1].plot(time, abs(Pressure - meanPress), label='Pressure',
                   **plot_param, color='midnightblue')
        ax[1].set_ylabel(ax_lab[2], fontsize=20)
        ax[1].set_xlabel(ax_lab[0], fontsize=20)
        ax[1].tick_params(labelcolor='black', labelsize=16, width=3)
    else:
        _ax.set_xlabel(ax_lab[0], fontsize=20)

    # other operations
    plt.autoscale(True,'both')
    plt.savefig('StressPlot.png')
    plt.show()
    plt.close(fig3d)
    print("[bold cyan]Pressure:[/]")
    print(" P = {:.1f} +/- {:.1f} kBar".format(Pressure[-lastConfig:].mean(),
                                               Pressure[-lastConfig:].std()
                                               )
          )


def _meanP(P):
    """
    Calculate the mean of the pressure from t=0 to a current time t

    Parameters
    ----------
    P : np.array
        the pressure data.

    Returns
    -------
    P_out : np.array
        mean of the pressure.

    """
    P_out = np.zeros(P.shape)
    for i in range(P.shape[0]):
        P_out[i] = P[0:i+1].mean()
    return P_out 

def plotPressCelldim(plotMeanStress=True, lastConfig=200, both=True,
                     repTime=False, tstep=5.0, plotDiff=False, dpi=150):
    """
    

    Parameters
    ----------
    plotMeanStress : TYPE, optional
        DESCRIPTION. The default is True.
    lastConfig : TYPE, optional
        DESCRIPTION. The default is 200.
    both : TYPE, optional
        DESCRIPTION. The default is True.
    repTime : bool, optional
        If True, fix time accumulation activation problem.
        The default is false.
    tstep : float, optional 
        The time step in ps.
    plotDiff : bool, optional
        If True abs(P - Pmean) will also be plotted on a new axe.
        The default is False.
        N.B.: To be only activated for when only interested in Stress (not 
            not in cell dim)
    Returns
    -------
    None.

    """
    time = []
    Stress = readSTRESS(time, name=r'STRESS')
    if repTime:
        time = repair_time(time=time, forStrs=True)
    time = time * tstep

    Pressure = np.array(CompPressure(Stress))

    if both:
        plotBoth(time, Pressure, plotMeanStress=plotMeanStress,
                 lastConfig=lastConfig, repTime=repTime, tstep=tstep, dpi=dpi)
    else:
        onlyPress(time, Pressure, plotMeanStress=plotMeanStress,
                 lastConfig=lastConfig, plotDiff=plotDiff, dpi=dpi)


#______________________________________________________________________________
#                                   For AEMD
#______________________________________________________________________________

linear = lambda x,a: a*x

def delayTime(t, att):
    """
     find the index of the inferior bound we want to set for the time.
     Will ulterly allow delaying the time from which the plot will begin at


    Parameters
    ----------
    t : np.array of float
        time  array.
    att : float
        Inferior bound of the time.

    Returns
    -------
    indice : the index 
        DESCRIPTION.

    """

    for ind, val in enumerate(t):
        if val >= att:
            indice = ind
            break
    return indice

def sliceTemp(start, end, Temp):
    """
    Set boundaries of the temperature for plotting purposes.

    Parameters
    ----------
    start : integer
        index of the inferior bound setted for the time in time array.
    end : TYPE
        index of the superior bound setted for the time in time array.
    Temp : np.array 
        Temperature.

    Returns
    -------
    np.array
       containting the temperature in this fixed range (of index

    """
    t = []
    for i in range(start, end):
        t.append(Temp[i])
    return np.array(t), len(t)

def kappa(tau, Lm, N_ions, cellVol):
    """
    Calculate the thermal conductivity for bulk system in AEMD

    Parameters
    ----------
    tau : float
        characteristic time for equilibrium
    Lm : float
        spatial periodicity in the direction of the heat flux.
    N_ions : integer
        The number of ions in the considered system (Can be different from
        the total of the global system when only regarding a portion of it)
    cellVol : float
        volume of the system.

    Returns
    -------
    k : float
        heat conductivity of the material.

    """
    k = (Lm**2 / (4 * np.pi**2)) * (3*N_ions * (kB / cellVol)) / tau 
    return k
     
def _meanT(T):
    T_out = np.zeros(T.shape)
    for i in range(T.shape[0]):
        T_out[i] = T[0:i+1].mean()
    return T_out

def plot_AEMD(t, T_1, T_2, filename, title, posi, N,
              Lm, cellVol, delay=False,
              wait=None, xlim=None, phase2=False, inplane=False,
              fontSizeLeg=15, plot_lw=2.0,
              setlimDT=False,
              limDT=(10,300),
              titlBlc1=r"T$_\mathrm{C}$", titleBlc2=r"T$_\mathrm{F}$",
              separate=False, svnameDeltaT:str="DeltaT.png"
              ):
    
    deltaT = T_1 - T_2
    
    # figure parameters
    fig_param = dict(layout='constrained',
                     edgecolor='white', facecolor='white',
                     dpi=150, linewidth=plot_lw)
    
    # settings for the plot
    if separate:
        fig1, ax1 = plt.subplots(figsize=(9, 3.75),**fig_param,sharex=True)
        ax = ax1,
    else:
        fig3d, ax = plt.subplots(figsize=(9, 7.5),nrows=2, ncols=1,
                                 **fig_param,sharex=True)
    
    # labels for the different plots of differents axis
    tup_plot = (titlBlc1, titleBlc2)
    
    # parameter for plots
    plot_param = dict(marker=None, linestyle=None,
                      linewidth=plot_lw, markersize=None)
    
    # parameters for the legend of the different axis
    dic_legend= dict(loc=0, shadow=False, fontsize=fontSizeLeg)

    #parameter for the labels of the different axis
    ax_lab = ((u'', r'T (K)'),
              (u'Time (ps)', r'$\Delta$T (K)'))
    
    t0 = 0
    if delay:
        t0 = wait

    ax[0].plot(t + t0 , T_1, label=tup_plot[0], **plot_param, color='tomato')
    ax[0].plot(t + t0, T_2, label=tup_plot[1], **plot_param)
    ax[0].set_ylabel(ax_lab[0][1], fontsize=20)
    #xticks
    ax[0].tick_params(labelcolor='black', labelsize=18, width=3,
                      axis='both', direction='in',
                   length=4.0)
    # setting of the legend
    ax[0].legend(**dic_legend)    
    ax[0].set_title(title)
    if xlim is not None:
        ax[0].set_xlim(xlim)
    else:
        plt.autoscale(True,'both')

    if separate:
        ax[0].set_xlabel(ax_lab[1][0], fontsize=20)
        # other operations
        plt.savefig(filename)
        plt.show()
        plt.close(fig1)
        
        fig2, ax2 = plt.subplots(figsize=(9, 3.75),**fig_param,sharex=True)
        ax = ax2, ax2
        

    # setting of the label to the axis
    ax[1].plot(t + t0, deltaT, label = r'$\Delta$T', **plot_param,
               color='orchid')
    ax[1].set_xlabel(ax_lab[1][0], fontsize=20)
    ax[1].set_ylabel(ax_lab[1][1], fontsize=20)
    ax[1].tick_params(labelcolor='black', labelsize=18, width=3,
                      axis='both', direction='in', length=4.0)
    
    if phase2:
        # fit expo decrois de dist
        popt, pcov = curve_fit(expoDecr, t, deltaT)
        ax[1].plot(t + t0, expoDecr(t, popt[0], popt[1]),
               label=r'fit $A \times \exp{(-b x)}$', color="#63b8ff", lw=1.5)
        ax[1].set_yscale('log')

        tau = 1/popt[1]
        dpopt = np.sqrt(pcov[1][1])
        dtau = dpopt / (popt[1] ** 2)
        print()
        if inplane:
            G = kappa(tau*1e-12, Lm, N, cellVol)
            dG = G * dtau / tau 

            msg = r'$\kappa$ = {0:3.2e} $\pm$ {1:3.2e}'
            msg = msg.format(G, dG) + ' W$\cdot$m$^{-1}$.K$^{-1}$'

        else:
            G = 3*N*kB/(4.0*tau*1e-12) * 1e9
            dG = G * dtau / tau
#            g = (G*1e-9)*Lm/cellVol * 1e-9
#            dg = dG/G * g
#            R = 1/G
#            dR = dG / (G**2)
            msg = (r'G = {:3.2e} $\pm$ {:3.2e}'.format(G, dG) +
                   r' nW$\cdot$K$^{-1}$')
#            msg2 = (r'R = {:.3e} $\pm$ {:.3e}'.format(R, dR) +
#                    r' K$\cdot$nW$^{-1}$')
#            msg3 =  (r'g = {:.3e} $\pm$ {:.3e}'.format(g, dg) +
#                     r' GW$\cdot$K$^{-1} \cdot$ m$^{-2}$')
#            ax[1].text(posi[0], posi[1], msg2,
#                       fontsize=16, color='k', transform=ax[1].transAxes
#                       )
#            ax[1].text(posi[0], posi[1] + 0.30, msg3,
#                       fontsize=16, color='k', transform=ax[1].transAxes)
        ax[1].legend()
        ax[1].text(posi[0], posi[1] + 0.10,
                   r'$\tau$ = {0:.2f} $\pm$ {1:.2f} ps'.format(tau, dtau),
                   fontsize=16, color='k', transform=ax[1].transAxes
                   )
        ax[1].text(posi[0], posi[1] + 0.20, msg,
                   fontsize=16, color='k', transform=ax[1].transAxes)


    
    if xlim is not None:
        ax[1].set_xlim(xlim)
    else:
        plt.autoscale(True,'both')
    if setlimDT:
        ax[1].set_ylim(limDT)
    if separate:
        ax[1].set_yscale('log')
        # other operations
        plt.savefig(svnameDeltaT)
        plt.show()
        plt.close(fig2)
        return None, None
    else:
        # other operations
        plt.savefig(filename)
        plt.show()
        plt.close(fig3d)
        return fig3d, ax

expoDecr =  lambda x, A, B : A*np.exp(-B*x)

def opTemp(time, T1, T2, timeLim=False, endTime=None, delay=False, wait=None):
    """
    Operation on temperature trajectory

    Parameters
    ----------
    time : np.array

    T1 : np.array
        Temperature in bloc 1.
    T2 : np.array
        Temperature in bloc 2.
    timeLim : bool, optional
        Low(time)-pass filter. The default is False.
    endTime : float, optional
        Cutting time associated with timeLim. The default is None.
    delay : bool, optional
        High(time)-pass filter. The default is False.
    wait : float, optional
        Cutting time associated with delay. The default is None.

    Returns
    -------
    timeOut : np.array
        Output (with reduced size).
    Tb1 : np.array
        Output bloc 1 temperature.
    Tb2 : np.array
        Output bloc 2 temperature.

    """
    if timeLim:
        timeOut = time[time <= endTime]
        sis = len(timeOut)
        index=0
        if delay:
            index = delayTime(timeOut, wait)
            timeOut = timeOut[timeOut >= wait] - timeOut[index]
            

        Tb1, _ = sliceTemp(index, sis, T1)
        Tb2, _ = sliceTemp(index, sis, T2)
        return timeOut, Tb1, Tb2
    else :
        return time, T1, T2
    
    
#______________________________________________________________________________
#                               DOS
#______________________________________________________________________________
def eDosToArray(line, data, fermi=None):
    """
    

    Parameters
    ----------
    line : str
        a line read from the file.
    data : list
        the list that will contain eigen values.
    fermi : int or None, optional
        if None then we are in the case of a vibrational analysis.
        The default is None.

    Returns
    -------
    None.

    """
    a = np.array(line.split()).astype(float)
    if fermi is not None: 
        data.append(a[1])
        if a[1] != a[-2]:
            data.append(a[-2])
    else:
         for val in a:
             if val > 0:
                 data.append(val)
                 
def gaussDos(x, mu, sigma, typ='e'):
    """
    gaussian function

    Parameters
    ----------
    x : float or int
        position at which the gaussian is evaluated.
    mu : float
        the mode of the gaussian.
    sigma : float
        half-heigh width.
    typ : str, optional
        'e' for electron and 'v' for vibration. The default is 'e'.

    Returns
    -------
    float

    """
    if typ.upper() == 'E':
        w = 1. / (sigma * np.sqrt(2*np.pi))
    elif typ.upper() == 'V':
        w=1
    return  w * np.exp(-1 * ((x - mu) / sigma)**2 )
    
def extract(filename, fileOut, prekeyWord, keyword, endword, skip, ncols,
            fermi, NbUnOcState=1
            
             ):
    """
    Extract eigenvalues (eDos or vDos) informations from the cpmd output file.

    Parameters
    ----------
    filename : str
        The name of the CPMD output file.
    fileOut : str
        The name of the file in which we want to copy() the important
        informations.
    prekeyWord : str
        Indicate from where (in output file) the interesting informations
        begin.
    keyword : str
        Indicate precisely the last (non empty) liine before the eigenvalues in
        the outupt file.
    endword : str
        word just after the eigenvalues.
    skip : bool
        If true then the line following keyword is ignored.
    ncols : int
        the number of columns in which the informations are written in CPMD
        output file.
    fermi : int or None
        if None then vibrational file else the considered file is for
        electronic eigenval.
    NbUnOcState : int, optional
        The number of unoccupied state energy we want to include our output
        data. The default is 2.

    Returns
    -------
    np.array
        the array of the eigenvalues.

    """
    data = []
    
    itrUnOcState=0
    ofi = open(filename, 'r')
    ofo = open(fileOut, 'w')

    while 1:
        line = ofi.readline().strip()
        
        if prekeyWord in line:
            while keyword not in line:
                line = ofi.readline().strip()
            if skip:
                ofi.readline()
            print(line)
            if fermi is not None:
                
                while endword not in line:
                    line = ofi.readline().strip().replace("NC", "0")
                    if '0.00000000' in line:
                        itrUnOcState += line.count('0.00000000')
                    if itrUnOcState <= NbUnOcState:
                        eDosToArray(line, data, fermi=fermi)
                        ofo.write(line+'\n')
                line = line.split()
                fermi = float(line[len(line)-2])
                print(f'[bold green]fermi[/] energy = {fermi:.2f} eV')
            else:
                
                while endword not in line:
                    line = ofi.readline().strip()
                    if endword in line:
                        break
                    else:
                        eDosToArray(line, data, fermi=fermi)
                        ofo.write(line+'\n')
            break
    ofi.close(); ofo.close
    
    return np.array(data).astype("float64"), fermi
    
#______________________________________________________________________________
#                            MSD Exciton
#______________________________________________________________________________

def generatorAxeIndex(ax, nc, nr):
    """
    generate axis
    Useful when dimensionality axis may vary (for example ax, ax[j] or ax[i,j])

    Parameters
    ----------
    ax : matplotlib.axes._axes.Axes
        axis created.
    nc : int
        number of columns in the figure.
    nr : int
        number of rows in the figure.

    Yields
    ------
    matplotlib.axes._axes.Axes
        Generate sequentially each axis contained in the figure.

    """
    if nc == 1 and nr == 1:
        yield ax
    elif max(nr, nc)>1 and min(nr, nc)==1:
        for i in range(max(nr, nc)):
            yield ax[i]
    else:
        for i in range(nr):
            for j in range(nc):
                yield ax[i, j]
                
#______________________________________________________________________________
#                           Transition moment (former)                        #
#______________________________________________________________________________
             
def parserTM0(path, dirname0='0', dirname1='1', dirname2='2'):
#This function has been depreciated the newer version is (parserTM)
    """
    Extract transition moment calculate different configurations

    Parameters
    ----------
    path : str
        Path to directories.
    dirname0 : str, optional
        Name of the directory corresponding to the Wavefunction (WF)
        optimization.
        The default is '0'.
    dirname1 : str, optional
        Name of the directory corresponding to the FEMD calculation of state
        occupation and energies.
        The default is '1'.
    dirname2 : str, optional
        Name of the directory corresponding to the calculation of 
        moment matrix elements.
        The default is '2'.

    Returns
    -------
    None.

    """

    dirn = []
    msg = '{}\t{}\t{}\n'
    
    # loop on the name of directories and files present in the path= './''
    for dirname, _, filenames in os.walk(path):
        
        # ignore dir that names end in '\\0' or '\\1' (for windows)
        # or '\0' and '\1' (for linux)
        if os.path.split(dirname) in [dirname0, dirname1] in dirname:
            continue
        else:
            # remaining dir names (in the present case) end by '\\2'
            if(dirname[-1] in dirname2 and dirname[-2] in '\\'):
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
 
    
#______________________________________________________________________________
#                           construct CPMD inputfile
#______________________________________________________________________________

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
#______________________________________________________________________________
#                   Transition moment calculation                             #
#______________________________________________________________________________
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
        e.g.:
            If the model input file is 'p3ht.in' for example,
            in the directory ./0/ (source), there is another directory '0' 
           (target) that contains the input file 'p3ht.in'.
    target : str, optional
        Name of the target directory (in source directory) that contains
        filename.
        The default is '0'.
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
    msg = '[bold green]{}[/] not found.'
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
                
                print("[bold green]%s[/] already exist !" %dest)
                print("Would you like to operate on it ([bold green]yes[/]" +
                " [bold purple]or[/] [bold red]no[/]) ?")
                overwrite = input(">> ")
                if overwrite.strip().upper() == 'NO':
                    continue
                else:
                    os.chdir(os.path.join(dest, target))
            elif compute.strip().upper() == "NONE":
                continue
        else:
            copytree(src=source, dst=dest)
            os.chdir(os.path.join(dest, target))
        print(">> The directory [bold green]%s[/] has been created." %dest[2:])
        
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
        dstFile.close()
        
    trajec.close()
    
def gauss(x, mu, height, sigma):
    return height * (1. / (sigma * np.sqrt(2*np.pi))) * np.exp(
                                                (-1./2.0) * ((x - mu) / 
                                                             sigma)**2 )

def Lorentz(x, mu, height, sigma):
    return height * (1. / (1 + ((x - mu) / (sigma / 2.) )**2 ))

def gaussBroadn(Eigen, minVal, maxVal, npoints, width, height=1.):
    """
    Perform gaussian broadening from EigenVal (eventually and oscillator
    strength) and return the spectrum data for plotting

    Parameters
    ----------
    Eigen : np.array
        np.array of eigen val.
    minVal : float
        minimum of the spectrum.
    maxVal : float

    npoints : int
        Desired number of points.
    height : float or np.array of float
        The height of the gaussian. the default is 1.0.
    width : float or np.array of float
        guassian width.

    Returns
    -------
    Enrg : np.array
        x-ax val of the spectrum.
    magnitude : np.array
        y-ax.

    """
    Enrg = np.linspace(minVal, maxVal, npoints)
    magnitude = np.sum(gauss(x=Enrg[:,np.newaxis], mu=Eigen, height=height,
                             sigma=width),
                       axis=1
                       )
    return Enrg, magnitude

def fmtMsg(line, linenumCols):
    out=[]
    for num in linenumCols:
        out.append(line[num])
    return tuple(out)
    

def parserTM(pathCalc='./', nameTargetDir='2', nameTargetFile="MATRIX.DAT",
             outFileName='transMoment.dat', numCols=(5, 6, 10)):
    path = os.getcwd()
    print('\n')
    print(' '*38, "[bold cyan]CWD:[/]")
    print("[bold yellow]%s[/]" %path)
    print('\n')
    
    dirn = []
    msg = '{}\t'*(len(numCols)-1)+'{}\n'
    nameTargetDir = '\\'+nameTargetDir
    
    # # loop on the name of directories and files present in the path= './''
    for dirname, _, filenames in os.walk(pathCalc):
        # ignore all dir that does not contains the targeted directory
        if(dirname.endswith(nameTargetDir)):
            dirn.append(dirname)
        else:
            continue
    
    # remove upperlevel dir that have same name as target dir
    unwanteDir = path.strip("\\")+nameTargetDir
    if unwanteDir in dirn:
        dirn.remove(unwanteDir)
      
    with open(outFileName, 'w') as file:
        for dirname in dirn:
            #os.path.join(path, dirname) contatenate path and dirname
            os.chdir(os.path.join(path, dirname))
            try:
                # first, test if MATRIX.DAT exist in current dir
                open(nameTargetFile, 'r')
            except OSError:
                # if not, close, print message & continue on followin dir
                print('No [bold green]{0:}[/] in [bold green]./{1:}[/]'.format(
                    nameTargetFile, dirname.split('\\')[-1])
                      )
                continue
            else:
                # if no exception (i.e MATRIX.DAT exist in current dir)
                with open(nameTargetFile, 'r') as matrix:
                    while 1:
                        line = matrix.readline()
                        if line: # True if line != '' (line not empty)
                            line2 = line
                        else:
                            matrix.close()
                            line2 = line2.split()
                            # write energies and trans mom in new file
                            file.write(msg.format(*fmtMsg(line2, numCols)))
                            break   
              
                
def pltSpctr(E, DeltaE, f12, oscStrength, w, minE, maxE,
             npoints, ax, pltColor, barwidth, barColor, lab):
        # gaussian broadening of the spectra
        Enrg, g=gaussBroadn(Eigen=E, minVal=minE, maxVal=maxE,
                            npoints=npoints, height=f12, width=w)
        
        ax.plot(Enrg, g / max(g), alpha=0.8, linewidth=2.0,
                color=pltColor,
                label='spctr %s' %lab)
        if oscStrength:  
            ax.bar(DeltaE /  q, f12 / max(f12), width=barwidth,
                   label=r"f$_{12}$ %s" %lab, color=barColor)              


def calcspctrDat(err, E_lumo, E_homo, tm, histNRJ):
    
    DeltaE = (E_lumo - E_homo + err) * q # J
    E = (E_lumo - E_homo + err)
    dip = tm * bohr**2
    
    #shift flexible 
    f12 = (2./3.) * (me / hbar**2) * (DeltaE - err*q) * dip
    if histNRJ:
        plt.figure(figsize=(6.5, 5.), dpi=200)
        plt.hist(E, density=False, bins=len(E))
        plt.xlabel("eV")
        plt.show()
    return E, DeltaE, f12, dip
                
def calcTau(E_homo, E_lumo, tm, D, T, err, minE, maxE, plot,
            pltColor, barColor, npoints, fontsize, barwidth, figsize, dpi,
            spltColor, sbarColor, oscStrength, superpose, sminE, smaxE,
            sE_homo, sE_lumo, stm, locLeg, xlim, ylim, labels, printLifeT,
            histNRJ, return_A21:bool=False
            ):
    
    

    print("\n"*2)
    print(" "*35+"[bold cyan]COMPUTING EXCITON LIFETIME[/]")
    
    D = D * 1e14 #cm2/s
    
    #err = Eg_exp - E_00

    E, DeltaE, f12, dip = calcspctrDat(err, E_lumo, E_homo, tm, histNRJ)
    
    if printLifeT:
        # A SECOND EXPRESSION OF A21
        A21 = 2 * (DeltaE**3) * (q**2) * dip / (3 * np.pi * hbar**4 *
                                                eps0 * c**3)
        np.savetxt("A21.dat", A21)
        meanA = np.mean(A21)
        #stdA = np.std(A21)
        tau = 1. / meanA
        #tau = np.mean(1. / A21)
        #meantau2 = np.mean(1/ sum(A21)**2)
        tau2 = tau**2
        meantau2 = 1. / np.mean((A21*A21))
        print((np.sqrt(abs(meantau2 - tau2))*1e12))
        #print(tau*tau)
        #Dtau = np.sqrt(meantau2 - tau*tau)
        #Dtau = tau * stdA/meanA
        # print("meanA {} stdA {} Dt(ps) {}".format(meanA, stdA, Dtau))
        # print("minA {} maxA {}".format(min(A21), max(A21)))
        #print(Dtau*1e12)
        print("[bold green]The exciton lifetime[/]")
        print("tau = %.2f [bold red]ps[/]" %(tau * 1e12))
        print("[bold green]The exciton diffusion length[/]")
        print('L = {:.2f} [bold red]nm[/]'.format(np.sqrt(D * tau)))
        print("\n"*2)
        print(" "*35+"[bold cyan]PLOT OF THE SPECTRUM[/]")
    
    if plot:
        w = convKToeV * T #
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figsize, dpi=dpi)
        
        pltSpctr(E=E, DeltaE=DeltaE, f12=f12, oscStrength=oscStrength, w=w,
                 minE=minE, maxE=maxE, npoints=npoints, ax=ax,
                 pltColor=pltColor, barwidth=barwidth, barColor=barColor,
                 lab=labels[0]
                 )
        
        # Second spectrum
        if superpose:
            sE, sDeltaE, sf12, _ = calcspctrDat(err, sE_lumo, sE_homo, stm,
                                                histNRJ)
            
            pltSpctr(sE, sDeltaE, sf12, oscStrength, w, sminE, smaxE,
                     npoints, ax, spltColor, barwidth, sbarColor,
                     lab=labels[1])
        
        ax.set_xlabel('$E$ (eV)')
        ax.set_ylabel(r'Intensity (arb. un.)')
        ax.tick_params(axis='both', direction='inout', length=4.0, width=1.0)
        
        print("The intensities are normalized regarding the "+
              "highest peak in the spectrum.")
        
        leg=dict(loc=0, shadow=False, fontsize=fontsize)
        
        if oscStrength or superpose:
            plt.legend(**leg)
        
            
        if xlim==None and ylim==None:
            plt.autoscale(True,'both')  
        elif not xlim==None:
            ax.set_xlim(xlim)
        elif not ylim==None:
            ax.set_ylim(ylim)
        else:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            
            
        plt.savefig('spectr.png', dpi=dpi)
        plt.show()
        plt.close(fig)
    if return_A21:
        return A21

#______________________________________________________________________________
#                               FixGeo
#______________________________________________________________________________
def distance(arr):
    """
    Calculate the norms of each column in an array 

    Parameters
    ----------
    arr : np.array
        The input array.

    Returns
    -------
    dist : np.array
        the norms.

    """
    dist = np.empty((arr.shape[1],))
    for i in range(arr.shape[1]):
        dist[i] = np.linalg.norm(arr[:, i])
    return dist


def transVec(cellVec1, cellVec2, cellVec3):
    """
    All possible translation (neighbours cells) of a cell given its
    lattice vectors

    Parameters
    ----------
    cellVec1 : np.array
        The first lattice vector.
    cellVec2 : np.array
        The second.
    cellVec3 : np.array
        The third.

    Returns
    -------
    trans_vec : np.array
        an array of all the possible translation vectors.
        array of shape (3, 28)

    """
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

def sliceArray(array:np.array, inf:float, sup:float, col:int=None):
    """
    Slice an array

    Parameters
    ----------
    array : np.array
        The array to be sliced.
    inf : float
        The smallest value.
    sup : float
        The highest value.
    col : int, optional
        The index of the column of reference. The default is None.

    Returns
    -------
    np.array

    """
    if col != None:
        msk1 = array[col, :] >= inf
        msk2 = array[col, :] <= sup
        msk = msk1 * msk2
        return array[:,msk]
    else:
        msk1 = array >= inf
        msk2 = array <= sup
        msk = msk1 * msk2
        return array[msk]
        
    
    
def finMaxNrj(array:np.array, inf:float, sup:float, target:int, col:int=None):
    """
    Find the max the coord of the max val in a given interval

    Parameters
    ----------
    array : np.array
        The data.
    col : int, optional
        The column in abscisse. The default is None.
    inf : float
        Smallest value to be considered in col.
    sup : float
        The highest.
    target : int
        The target column.

    Returns
    -------
    TYPE
        DESCRIPTION.
    maxtarg : TYPE
        DESCRIPTION.

    """
    array = sliceArray(array, inf, sup, col)
    maxtarg = max(array[target, :])
    indMax = list(array[target, :]).index(maxtarg)
    return array[col, indMax], maxtarg

def cartToSphr(x, y, z):
    """
    Transform cartesian coordinates to spherical coordinates

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.

    Returns
    -------
    theta : TYPE
        DESCRIPTION.
    phi : TYPE
        DESCRIPTION.

    """
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

def norm(arr:np.array, arr2:np.array=None):
    """
    Calculate the distance between two vector
    
    If only one vector is inputted, this function calculate its norm.
    Parameters
    ----------
    arr : np.array
        the first array.
    arr2 : np.array, optional
        The second array.
        The default is None.
        
    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    if (arr2 == None):
        return np.sqrt(sum(arr*arr))
    else:
        return np.sqrt(abs(sum(arr*arr2)))

def sphericalVec(theta:float, phi:float):
    """
    Calculate the basis vectors of the spherical coord

    Parameters
    ----------
    theta : float
        polar angle.
    phi : float
        azimutal angle.

    Returns
    -------
    e_r : TYPE
        radial unit vector.
    e_theta : TYPE
        unit vector along theta.
    e_phi : TYPE
        unit vector along phi.

    """
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

def reorganize_data(*data, min_val:float=0.0, maxval:float=0.0):
    arr = sliceArray(data[0],inf=min_val, sup=maxval)
    maxv = arr.max()
    idx_max = list(data[0]).index(maxv)
    new_dat = []
    for i in range(idx_max+1, data[1].shape[0]):
        new_dat.append(data[1][i])
    for j in range(idx_max+1):
        new_dat.append(data[1][j])
    return np.array(new_dat)
