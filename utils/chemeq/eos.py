"""
Utility programs for computing and plotting results of the equation of state,
specifically the abundances of 92 neutral and ionized atoms and 503 molecular species,
for s specified temperature and one state parameter (in principle total particle density, 
mass density, or electron density).
At the moment, there are two limitations, which can be easily avoided, but on the
cost of having the codes and corresponding input heavier.
- only N_tot is accepted as a state parameter
- solar abundances are assumed (after Asplund 2009)

There are four individual routines:

1) state - calls the F90 program chemeq that performs the actual calculations,
(essentially using several slightly modified subroutines from Synspec), 
and makes two simple  plots, one for atoms, one for molecules (by calling routines
atoms and molec, which call the basic plotting routine plabund.
The Fortran program chemeq produces the standard output (similar to the one in Synspec)
and two files fort.78 (for atoms/ions) and fort.79 (for molecules), called EOS files.

2) plabund - plots and annotates the abundances of either atoms/ions, or molecules.

3) atoms - just a shortcut for calling plabund for atomic species with a default setup

4) molec - the same for molecular species.

Examples:
--------

1) To solve the atomic/molecular equilibrium for given parameters (T,Ntot),
and to produce initial plots:

>>> import eos
>>> eos.state(t=3000,ntot=1.e13)

2) To replot molecular abundances with a different marker shape and size:
 
>>> x,y,lab = eos.plabund(mark='x',msize=4)

This will plot a rather cluttered graph, so one will need to enlarge it and
select a part of interest. 

3) to overplot the current plot with previously computed file fort.79, 
called "prev.79", and without additional annotation:

>>> x,y,lab = eos.plabund(file='prev.79',iover=1,annot=0)   

"""


################################################################################

import matplotlib.pyplot as plt
import numpy as np
import os


################################################################################


def state(t=4000,ntot=1.e12):
    
    """
    Parameters:
    - t - temperatuture (in K)
    - ntot - total particle numbere density [cm^-3]
    """

    f = open('tmp.5','w')
    f.write('%f %e \n' %(t,ntot))
    f.close()

    a = "./chemeq < tmp.5 > run.log"
#
#   This is to be replaced by the actual location and name of the
#   use executable 
    print(a)
    os.system(a)

    molec()
    plt.xlim(0,100)
    lne = np.log10(ntot)
    titl = 'T = ' + str(t) + 'K,  log N_tot = ' + str(lne)
    plt.title(titl)

    atoms()
    plt.title(titl)


################################################################################


def plabund(file='fort.79',ix=0,iy=3,mark='o',msize=2,iover=0,annot=1):

    """
    Plots and annotates the abundances of either atoms/ions, or molecules.
    the plots display logarithm of abundance as a function of atomic number (for atoms/ions),
    or as a function of Tsuji index (for molecules). The routine recognizes which is which,
    because the input file for  atoms/ions has 6 columns, while that for molecules just 4 columns.
    To avoid clutter, the molecular plots is limited to the first 100 molecular species, but
    this can be changed by setting plt.xlim appropriately.
    The routine, called stand-alone, does not provide a title of the plot (unlike when
    it is called from routine state.

    Parameters:
    file: str, optional 
         input file (result of a run of chemeq)
         default: "fort.79' (last computed molecular EOS file
    ix:  integer, optional
         column for the x-coordinate of the plos
         default: 0
    iy:  integer, optional
         column of the y-coordinate
    mark: str, optional
         marker for the abundance plot
         default: 'o'
    msize: integer, optional
         size of the marker
         default: 2
    iover: integer, optional
         if 0, new plot
         else: overplot
    annot, integer, optional
         if 1, markers are annotated by a chemical lable
         else: no annotationo

    Returns; X,y,lab - the arrays of x, y, coordinates, and an
    array of chemical labels.
    """


    f = open(file)
    lines = f.readlines()

    i=0
    x = []
    y = []
    lab=[]
    while i < len(lines):
        xlin=lines[i]
        xlin = xlin.strip()
        val  = xlin.split()
        x.append(float(val[ix]))
        y.append(float(val[iy]))
        lab.append(val[1])
        i    += 1

    n = i-1
    if iover == 0:  plt.figure()
    plt.plot(x,y,marker=mark,ls='',ms=msize)
    if iy == 3: plt.xlabel('Tsuji index')
    if iy >  3: plt.xlabel('atomic number')
    plt.ylabel('log(abundance)')
   
    if annot == 1:
        if iy < 5:
            for i in range(0,n): plt.annotate(lab[i],xy=(x[i],y[i]),ha='center',size='smalll',rotation=30)
        if iy == 5: 
            for i in range(0,n): plt.annotate(lab[i]+'+',xy=(x[i],y[i]),ha='center',size='x-small',rotation=30,color='r')
    
    return(x,y,lab)

################################################################################

def atoms(file='fort.78',iover=0):

    x,y,lab = plabund(file,iy=4)
    a,y,lab = plabund(file,iy=5,iover=1)

################################################################################

def molec(file='fort.79',iover=0):

    x,y,lab = plabund(file,iover=0)


################################################################################



