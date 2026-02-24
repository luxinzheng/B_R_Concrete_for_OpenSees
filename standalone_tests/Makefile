
FC = gfortran
FFLAGS ?= -O2 -std=f2008 -Wall -Wextra -fimplicit-none -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
LDFLAGS ?=

OBJS = forumat.o br_test_params.o br_test_utils.o test_runner.o

all: test_runner

forumat.o: ../forumat.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

test_runner: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $@ $(LDFLAGS)

run: test_runner
	./test_runner

clean:
	rm -f *.o *.mod test_runner *.csv
