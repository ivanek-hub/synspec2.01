The directory contains input data for computing a simple test opacity table.
The table has 10,000 wavelength points between 900 and 11,000 A,
for 3 temperatures between 3 and 15 kK (spaced logarithmically),
and 3 densities, between 10<sup>-12</sup> and 10<sup>-6</sup> g cm<sup>-3</sup> 
(again spaced logarithmically).
The resulting opacity table is in the file *optes.dat*, as stipulated 
in the input file *fort.2*.

The input data assume that the binary line lists are being used, so that
they must exist. In the case the binary lists are not used, the user
should make links to the ASCII files, and to make corresponding changes in
the file *fort.3*, namely

'gfATO.dat'          0 <br>
'gfMOL.dat'          0    8000. <br>
'gfTiO.dat'          0    5000.

in which case the text (ASCII) format of line lists is used. 

For actual applications, more extend wavelength range, and a higher resolutions in temperature and
density are obviously needed. To show a construction of such opacity table (with better, but still
quite modest resolutions), we include two optional input files, to be copied to *fort.2*
(analogous to the test cases considered in Paper IV):

*optab11_7f.2* -- opacity table with all chemical elements, with higher resolutions; <br>
*optab11_7hmc.2* -- opacity table with opacity of H, Mg, and Ca removed.
