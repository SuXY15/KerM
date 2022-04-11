PROGRAM := mixing_test
FC      := gfortran
FCFLAGS := -g # -O
# FCFLAGS += -ffree-line-length-none
# FCFLAGS += -g -fbacktrace -fno-align-commons -fbounds-check
FLFLAGS := 

SRCS_f   := $(wildcard *.f)
OBJS_f   := $(SRCS_f:%.f=%.o)
SRCS_f90 := $(wildcard *.f90)
OBJS_f90 := $(SRCS_f90:%.f90=%.o)

default:
	make $(PROGRAM)
	rm -f *.o *.mod *.out

$(PROGRAM):$(OBJS_f) $(OBJS_f90)
	$(FC) $(FLFLAGS) -o $@ $^

$(OBJS_f):%.o:%.f
	$(FC) -c $(FCFLAGS) -o $@ $<

$(OBJS_f90):%.o:%.f90
	$(FC) -c $(FCFLAGS) -o $@ $<

#------------  targets  ----------------------------------------------
mixing_test.o: emst_subs.o
emst.o: emst_subs.o

clean:
	rm -f *.o *.mod *.out