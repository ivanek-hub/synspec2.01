#
#	makefile for synspec
#

F90 = gfortran 
OPT = -Ofast -fno-automatic #-ftree-vectorize

SOURCE = synspec2.01.f 

SOURCE_MODULES = syn_modules.f syn_alloc.f

MODS = accura.mod  heprf.mod   modelp.mod  optabl.mod  synthp.mod \
       allarn.mod  hydprf.mod  molist.mod  params.mod  topdat.mod \
       eospar.mod  lindat.mod  opadat.mod  photcs.mod  wincom.mod

OBJS = $(SOURCE_MODULES:.f=.o) 

all: synspec2.01


# 2nd and final step - compiling and linking
synspec2.01:  ${OBJS} ${SOURCE}
	${F90}  -o synspec2.01 ${OPT}  ${SOURCE}

#1st step - compiling the modules to get .o and .mod files
%.o : %.f 
	${F90} -c ${OPT} ${SOURCE_MODULES}

clean: 
	rm -f synspec2.01 *.o *.mod


