### 1. EXO2SYN -- EXOMOL to SYNSPEC transformation of data

Contains the following programs, described in detail in the Manual:

- list.f: the program extracts line data for particular isotopologue of a molecule from
the EXOMOL database, and outputs them in the SYNSPC format.

- reverse.f: a simple program that reverses the order of wavelengths from decreasing
(in EXOMOL) to increasing (in Synspec)

- merge.f: merging two line lists

- subdirectory *pfsum_data*, containing a collection of input data for pfsum.

A collection of scripts and input data  for building the EXOMOL-based line lists using the above programs can be found [here](https://github.com/callendeprieto/exomol).

 