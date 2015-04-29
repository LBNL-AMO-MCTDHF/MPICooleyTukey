
include DFFTPACK/Dfftfiles.Inc

include Makefile.header


default: ft_test mpi_test 

clean:
	rm *.o *.a *genmod* *.mod *~ mpi_test ft_test; cd DFFTPACK; make clean

TESTSRC = cooleytukey.o ft_test.o ftcore.o cooleytukey_shared.o

MPISRC = cooleytukey_mpi.o mpi_test.o ftcore.o MPI.o cooleytukey_shared.o

ft_test:  $(TESTSRC) $(DFFTCHECK)
	$(F90) $(LOADFLAGS) -o ft_test $(TESTSRC) $(LIBRARIES)

mpi_test:  $(MPISRC) $(DFFTCHECK)
	$(F90) $(LOADFLAGS) -o mpi_test $(MPISRC) $(LIBRARIES)

cooleytukey.o:  cooleytukey.f90
	$(F90) $(FFLAGS)   -c cooleytukey.f90 -I$(MYINCLUDE) 

cooleytukey_shared.o:  cooleytukey_shared.f90
	$(F90) $(FFLAGS)   -c cooleytukey_shared.f90 -I$(MYINCLUDE) 

cooleytukey_mpi.o:  cooleytukey_mpi.f90
	$(F90) $(FFLAGS)   -c cooleytukey_mpi.f90 -I$(MYINCLUDE) 

MPI.o:  MPI.f90
	$(F90) $(FFLAGS) -c MPI.f90 -I$(MYINCLUDE)


ftcore.o:  ftcore.F90
	$(F90) $(FFLAGS) $(FFTFLAG)  -c ftcore.F90 -I$(MYINCLUDE) 


ft_test.o:  ft_test.f90
	$(F90) $(FFLAGS)   -c ft_test.f90 -I$(MYINCLUDE) 

mpi_test.o:  mpi_test.f90
	$(F90) $(FFLAGS)   -c mpi_test.f90 -I$(MYINCLUDE) 

DFFTPACK/passb.o:
	cd DFFTPACK; make



