FC = gfortran
FCFLAGS = -g -fcheck=all -Warray-bounds

# source files and objects
SRCS = donnees.f90 sequentiel.f90 projet.f90
# program name
PROGRAM = run

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

project.o:sequentiel.o
sequentiel.o:donnees.o

clean:
	rm -f *.o *.mod run
