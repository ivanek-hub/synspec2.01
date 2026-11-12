# SYNSPEC2.01

## Description

The depository contains the spectrum synthesis program SYNSPEC, version 2.01. It allows the calculation of synthetic spectra of stars, planets, and disks, from model atmospheres.The code has been migrated to modern FORTRAN, and supplemented with a suite of utilities that allow users to build compatible molecular line lists from EXOMOL data.


At the coding level, previous versions of the code, written in FORTRAN77, 
faced a number of problems on modern computers.
The program was  therefore converted to modern FORTRAN. 
The upgrade enables a more efficient handling and bookkeeping of large arrays needed 
to store atomic, and in particular molecular, line data needed for generating
synthetic spectra of cool objects.


Besides the SYNSPEC itself, we also present a number of utility programs, such as the 
suite of routines to generate molecular line lists in the SYNSPEC format from EXOMOL 
data, and a short program CHEMEQ for a fast exploration of a solution
of chemical equilibrium under various  conditions, and with various molecular
data. This program is also supplemented by a Python program to display the
results in a graphical way.



## Documentation

The basic document for the present version of SYNSPEC is a paper submitted to Astronomy and Astrophysics (Hubeny, Allende Prieto, and Osorio, 2027 - refereed to as Paper V) that contains a description of basic physics. It is accompanied by an additional document, called Manual, that contains a technical description and user's guide. Both documents are present in the directory *docs*.

For convenience, the directory *docs* also presents a four-volume series of papers which describe previous versions of SYNSPEC, as well as stellar atmosphere code TUSTY. Papers I–III cover TLUSTY 205 and SYNSPEC 51. Paper IV (Hubeny et al. 2021, arXiv:2104.02829) covers the last FORTRAN77 version of SYNSPEC (54). Although the most recent paper together with the Manual provide a comprehensive description of the current SYNSPEC, the previous papers still contain useful information about the basic physics, atomic data files, etc.


## Repository layout

There are five basic directories:

* data <br> 
Contains a collection of needed atomic data files, partition functions for atoms and molecules, and molecular equilibrium constants

* docs <br>
Contains the documentation, as specified above

* source <br>
Contains the source code, in three separate files: <br>
*synspec2.01.f* <br>
*syn_modules.f* <br>
*syn_alloc.f* <br>
together with the *makefile* file

* tests <br>
Contains four different test cases that illustrate a usage of Synspec in various applications. For each test case we present the script to run the test, the input files, and the subdirectory with the results of the test that can be compared to the result the user obtains, to verify a proper installation of the program and the data.


* utils <br>
Contains three groups of utility programs, described inside the directory.
