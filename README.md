> # SYNSPEC2.01

## Description

This repository contains the spectrum synthesis program SYNSPEC, version 2.01. It allows the calculation of synthetic spectra of stars, planets, and disks, from model atmospheres.The code has been migrated to modern FORTRAN, and supplemented with a suite of utilities that allow users to build compatible molecular line lists from EXOMOL data.


At the coding level, previous versions of the code, written in FORTRAN77, 
faced a number of problems on modern computers.
The program was  therefore converted to modern FORTRAN. 
The upgrade enables a more efficient handling and bookkeeping of large arrays needed 
to store atomic, and in particular molecular, line data needed for generating
synthetic spectra of cool objects.


Besides the SYNSPEC, we also present a number of utility programs, such as the 
suite of routines to generate molecular line lists in the SYNSPEC format from EXOMOL 
data, and a short program CHEMEQ for a fast exploration of a solution
of chemical equilibrium under various  conditions, and with various molecular
data. This program is also supplemented by a Python program to display the
results in a graphical way.



## Documentation

The basic document for the present version of SYNSPEC is a paper submitted to Astronomy and Astrophysics (Hubeny, Allende Prieto, and Osorio, 2027 - referred to as Paper V) that contains a description of the basic physics. It is accompanied by an additional document, called the Manual, that contains a technical description and the user's guide. Both documents are present in the directory *docs*.

For convenience, the directory *docs* also presents a four-volume series of papers which describe previous versions of SYNSPEC, and the stellar atmosphere code TUSTY. Papers I–III cover TLUSTY 205 and SYNSPEC 51. Paper IV (Hubeny et al. 2021, arXiv:2104.02829) covers the last FORTRAN77 version of SYNSPEC (54). Although the most recent paper together with the Manual provide a comprehensive description of the current SYNSPEC, the previous papers still contain useful information about the basic physics, atomic data files, etc.


## Repository layout

There are five basic directories:

* data <br> 
Contains a collection of required atomic data files, partition functions for atoms and molecules, and molecular equilibrium constants

* docs <br>
Contains the documentation, as specified above

* source <br>
Contains the source code, in three separate files, plus the Make file: <br>
> *synspec2.01.f* - main part: a collection of executive routines <br>
> *syn_modules.f* - a collection of modules <br>
> *syn_alloc.f* - a collection of array allocation routines <br>
> *makefile* 

* tests <br>
Contains four different test cases that illustrate the usage of Synspec in various applications. For each test case we present a script to run the test, the input files, and the subdirectory with the results of the test that can be compared to the result the user's own results, to verify the proper installation of the program and the data. <br>
> bstar - spectrum for a NLTE metal line-blanketed model atmosphere <br>
> hhe - spectrum for  a simple H-He NLTE model atmosphere <br>
> molec - spectrum of a cool Kurucz model atmosphere with multiple molecular line lists <br>
> optab - construction of a simple opacity table


* utils <br>
Contains three groups of utility programs in individual subdirectories: <br>
> chemeq - an exploratory equation-of-state solution, plus a Python program for the corresponding graphical display <br>
> exo2syn - a suite of programs to generate isotope-averaged molecular partition functions and line lists  from EXOMOL data <br>
> listhandl - short programs for handling (reducing) line lists


## Line lists

Since line lists are generally rather large files, they cannot be accessed directly from this repository, but instead from another source.

The test cases described here in directories *tests/hhe* and *tests/bstar* refer to a short atomic line list called *linelist_test.19*, which is a part of directory *data* here. But for realistic calculations of synthetic spectra, one needs a much more extended atomic line list.
Similarly, the test case described in *test/optab* uses this atomic line list, and the standard molecular line lists (containing all important molecules except TiO), plus a special list for TiO. 
These two, plus the line list for TiO, are available from 
[here](https://cloud.iac.es/index.php/s/6x5mRhXj7Ekyc9H?dir=/linelists/general/2022).


The line lists for individual molecules, used in the test case *tests/molec* are are available
[here](https://cloud.iac.es/index.php/s/6x5mRhXj7Ekyc9H?dir=/mollist).

Note that the files in these sites have somewhat different names then those referred to in the test cases, so the file *fort.3* in the test cases *tests/optab* should be modified, or the downloaded files renamed, accordingly.



 ## Graphical user interface


* *SYNPLOT* -- is an IDL wrapper
to prepare the necessary input files, run the code, possibly perform rotational
and instrumental convolutions, plot the results, and optionally provides
graphical line identifications. <br>
The  new variant of [SYNPLOT](https://github.com/argosun/synplot)
works in IDL and  GDL.



* *SYNPLE* -- is a Python wrapper. It is similar to *SYNPLOT*, with several important differences.
[SYNPLE](https://github.com/callendeprieto/synple)
 can take as input *MARCS*, Kurucz, *PHOENIX*, and *TLUSTY* model atmospheres, and is designed to be useful to researchers that are not experts in stellar atmospheres and radiative transfer. It facilitates  performing massive spectral synthesis calculations for large grids, parallelizing the calculations over multiple cores, using the Python libraries *multiprocessing* or *ray*,  and multiple processors, setting up scripts to run SYNSPEC  on distributed-memory  computers.
It includes as well a routine for fitting observed spectra using precomputed grids of
synthetic spectra. 
