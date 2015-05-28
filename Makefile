#
PREFIX = ./
LIBDIR = $(PREFIX)/lib
BINDIR = $(PREFIX)/bin
#
#######         For Linux              #########
FPATH      = $(shell which gfortran)
########        Compiling source       ##########

forcas:
	#f2py -m forcas -h forcas.pyf forcas.f --overwrite-signature 
	#f2py --f90exec=$(FPATH) -c forcas.pyf forcas.f
	rm -f *.o
	rm -f forcas.so
	f2py -c forcas.f90 -m forcas



install: installdirs
	install -m 0644 forcas.so $(LIBDIR)


installdirs:
	mkdir -p $(LIBDIR)


#################################################
########            Cleaning             ########
#################################################
clean:
	rm -f *.o
	rm -f forcas.pyf
	rm -f forcas.so
#
