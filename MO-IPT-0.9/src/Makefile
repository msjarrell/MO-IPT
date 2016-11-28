#fort = mpif90 -check all
fort = mpif90  
RM = rm -f
OBJ =  main.o  funct.o makegrid.o ksum.o 
MKLPATH= /opt/intel/mkl/lib/intel64/ 
MKL_LIBS=-L$(MKLPATH) -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5 -lpthread_nonshared -lfftw3 $(FFLAGS) $(SPRNGLIB) $(LFLAGS)

../MO_IPT.run : $(OBJ)
	$(fort) -O3 -o ../MO_IPT.run $(OBJ) $(MKL_LIBS)

main.o : main.f90
	$(fort) -O3 -c -o main.o main.f90

makegrid.o : makegrid.f90 
	$(fort) -O3 -c -o makegrid.o makegrid.f90 

funct.o : funct.f 
	$(fort) -O3 -c -o funct.o funct.f 

ksum.o : ksum.f90
	$(fort) -O3 -c -o ksum.o ksum.f90

clean:
	$(RM) $(OBJ) ../MO_IPT.run 
