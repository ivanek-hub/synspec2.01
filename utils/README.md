
## Utility programs

There are three groups of utility programs, each written in FORTRAN 90/95:

### 1. EXO2SYN -- EXOMOL to SYNSPEC transformation of data

Contains the following programs, described in the Manual):

- list.f: the program extracts line data for particular isotopologue of a molecule from
the EXOMOL database, and outputs them in the SYNSPC format.

- reverse.f: a simple program that reverses the order of wavelengths from decreasing
(in EXOMOL) to increasing (in Synspec)

- merge.f: merging two line lists

- subdirectory pfsum_data, containing a collection of input data for pfsum.

### 2. EOSPROG -- equation of state programs and data

Contains the following programs, described in Paper V:

- chemeq.f: chemical eqilibrium, or equation of state, solution
- eos_modules.f: auxiliary program for chemeg
- eos_alloc.f:   auxiliary program for chemeq
- isotops:       file with solar isotope fractions, needed for chemeq.f

- eos.py: a Python program for displaying the results of chemeq.f. User instructions
and explanations are given at the header of the program.


### 3. LISTHANDL -- handling of line lists

Contains 3 simple line-list-handling programs, described in Paper IV:

- redulis:  reduction of an original line list by removing weak lines
- list2bin: transformation of a line list to a binary format
- list2mini: transformation of a molecular line list to a "minimalistic" format (see Paper V).
