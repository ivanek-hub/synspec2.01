The directory contains input data for the test case,
computed using the included script file *R1*.
It demonstrates working with multiple molecular line lists, referred to in the file *fort.3*.
The model atmosphere is a  Kurucz model with
Teff=3500 K, log g = 4, for a wavelength range
between 3000 and 9000 A.

This example also shows the case where the previously used input file
*fort.55* is missing, and all keyword parameters are contained in the file 
*kwords*. Unlike previous versions of Synspec, no keyword parameter that
would indicate that the input model atmosphere is a Kurucz one
is needed, since Synspec now recognizes it automatically.

The input file *fort.3* invoked many individual molecular line lists. These are not a part of the standard distribution since many files are too big to be transferred through GitHub. They can be downloaded separately, and be placed to a special directory *$BASE/mollists*.
Analogously. the line list *gfATO.dat* (the full atomic line list), and *gtTiO.bin* should be located in a special directory *$BASE/linelist*, as in the previous release of Synspec.

The run also generates output files fort.12 (list of included atomic lines),
and fort.15 (the same for molecular lines). But they are not included
in the directory *results*, since they are quite big.
