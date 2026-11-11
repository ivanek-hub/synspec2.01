The directory contains input files for a test case of computing
a very short chunk of synthetic spectrum (1400-1410 A) 
for a NLTE metal line-blanketed model atmosphere from 
the BSTAR2007 grid (Lanz & Hubeny 2007), a model with 20,000 K, and log *g* = 4.

One can use the included script file *R1*.  
Notice that the command 

*cp fort.55.lin fort.55* 

is commented out, so the file *fort.55* is absent. This demonstrates the new system of input the keyword parameter, which are now
taken from the file *kwords* or, if not specified there, assume their default values.
For details, refer to the Manual. To recover the old system, one uncomments the line

*cp fort.55.lin fort.55*

and delete the file *kwords*.

The results of both runs are very similar, but small differences appear
due to different settings of the characteristic depth
(depth IDSTD=34 in the case of using *fort.55*, and default
values of 3 characteristic depth if IDSTD has its default
value IDSTD=0 -- see the Manual.

This test case also demonstrate the using of additional input file *dimens*, for resetting the dimensions of some arrays related to line lists. Finally, since the inout model is a NLTE metal line-blanketed model atmosphere, one needs an input file *superlevels* that determines the association of teh actual levels to the "superlevels"for which the NLTE level populations are given.
