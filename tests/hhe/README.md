The directory contains input files for generating a sample synthetic
spectrum of the model hhe35nl, i.e. a simple H-He NLTE model
atmosphere for *T*<sub>eff</sub>= 35,000 K, log *g* = 4, that was computed
as one of the test cases for TLUSTY, (see Paper IV).

One can use the included script file *R1*, or one can delete
or rename file fort.55, so the keyword parameters are taken
only from file *kwords*.

One can also use file *dimens*, e.g. the one from the test case *bstar*.
Without using it, the results are the same, but the run takes more memory.