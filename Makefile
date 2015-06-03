#
PREFIX = ./
FPATH      = $(shell which gfortran)
########        Compiling source       ##########

forcas:

	rm -f *.o
	rm -f forcas.so

	f2py --f90exec=$(FPATH) -c forcas.f90 -m forcas	
#	f2py -c forcas.f90 -m forcas


#################################################
########            Cleaning             ########
#################################################
clean:
	rm -f *.o
	rm -f forcas.pyf
	rm -f forcas.so
#
