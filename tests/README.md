This directory is composed of several subdirectories that contain scripts
to generate various typical test cases. The scripts assume that
prior to running the tests, the user sets the environment variable
BASE to point out to the parent directory that contains
subdirectory _./src_ with the source code and its executable, and
subdirectory _./data_ with a collection of atomic and molecular
data. Both subdirectories are a part of the standard distribution.

The _./data_ directory contains atomic data for several atoms/ions as
described in the Manual, data for atomic and molecular partition functions,
nd equilibrium constants for molecules. But the GitHub repositories DO NOT
contain atomic and molecular line lists, which are too big to be distributed
through GitHub. Some test cases (**hhe**, **bstar**) use a very short atomic line list
that is a part of the _./data_, but for a general spectrum synthesis the user should
download relevant line lists from different sources.

For example, if the parent directory is called _~/synspec_,
and contain subdirectories _./src_, and _./data_ extracted from the
corresponding repositories of the GitHub site, the environment
variable BASE is set (in zsh) as

_export BASE=~/synspec_