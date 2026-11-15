### EXO2SYN -- 

Contains the following programs which serve to transfer EXOMOL data (partition functions and line data, including averaging over isotopologues) to  the SYNSPEC-compatible  format. They are described in detail in the Manual:

- *pfsum.f*: program for averaging EXOMOL partition functions over isotopologues, assuming the solar isotope composition.

- *list.f*: the program extracts line data for particular isotopologue of a molecule from
the EXOMOL database, and outputs them in the SYNSPC format.

- *merge.f*: merging two line lists, used in this context for merging the line data for individual isotopologues of the given molecule.

- *reverse.f*: a simple program that reverses the order of wavelengths from decreasing
(in EXOMOL) to increasing (in Synspec)

In addition to the programs, this repository also contains:

- file *isotops*, specifying solar isotope composition, after Lodders (2020).
 
- subdirectory *pfsum_data*, containing a collection of input data for pfsum.

- A collection of scripts and input data  for building the EXOMOL-based line lists using the above programs, that can be found [here](https://github.com/callendeprieto/exomol).

 
