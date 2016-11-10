PROG =	run

SRCS =	donnees.f90 mod_preprocess.f90 mod_physique.f90 mod_fonction.f90 main.f90

OBJS =	donnees.o mod_preprocess.o mod_physique.o mod_fonction.o main.o

LIBS =	

CC = cc
CFLAGS = -O
FC = f77
FFLAGS = -O 
#F90 = ifort
F90 = gfortran
F90FLAGS = -O3 -march=native 
#F90FLAGS = -O3 -march=native -ftz -parallel 
#F90FLAGS = -O3 -march=native -ffpe-trap=invalid
#F90FLAGS = -O0 -pedantic -Wall -ffpe-trap=invalid,zero,overflow,underflow -g -fcheck=all -fbacktrace -fdump-core
#F90FLAGS = -g -O0 -traceback -check all -ftrapuv -fpe0 -qopenmp
#LDFLAGS = -parallel 
LDFLAGS =

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
	make cleanmodres

clean:
	rm -f $(PROG) $(OBJS) *.mod *~ fichier/*.dat gmon.out analysis.txt test_gprof

cleanmodres:
	rm -f $(OBJS) *.mod *~ fichier/*.dat 
	
gprof: 
	$(F90) -pg -O0 donnees.f90 mod_fonction.f90 mod_physique.f90 mod_preprocess.f90 main.f90 -o test_gprof 
	make cleanmodres
	./test_gprof
	gprof -l test_gprof gmon.out > analysis.txt

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

mod_fonction.o: donnees.o
mod_preprocess.o: donnees.o mod_fonction.o
mod_physique.o: donnees.o mod_fonction.o
main.o: mod_fonction.o mod_physique.o mod_preprocess.o donnees.o
