### CHEMEQ  

The program equation for solving the chemical equilibrium equations (the equation of state), together with its auxiliary programs. They are descibed in detail in Paper V.

- *chemeq.f*: chemical eqilibrium (equation of state) solver

- *eos_modules.f*: auxiliary program for *chemeg* - a shortened file *syn_modules* of the *src* directory - collection of modules

- *eos_alloc.f*:   auxiliary program for *chemeq* - a collection of allocation routine

- *isotops*:       file with solar isotope fractions, needed for *chemeq*

- *eos.py*: a Python program for displaying the results of *chemeq*. User instructions
and explanations are given at the header of the program.
