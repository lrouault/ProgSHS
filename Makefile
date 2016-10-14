FC = gfortran
FCFLAGS = -g -fcheck=all -Warray-bounds

# source files and objects
SRCS = sequentiel.f90 projet.f90
# program name
PROGRAM = test

all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

project.o:sequentiel.o

clean:
	rm -f *.o *.mod test
