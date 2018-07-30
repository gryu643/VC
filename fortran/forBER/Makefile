# Makefile for colMat.f90 noiseExist

CC = gfortran
PROGRAM = CALmod PPLmod usePPLmod
CFLAGS = -O3
OBJS = CALmod.o PPLmod.o usePPLmod.o
MOD = pplmod.mod calmod.mod
LAPACK = -llapack

usePPLmod: CALmod.o PPLmod.o usePPLmod.o
	$(CC) $^ -o $@ $(LAPACK)

vc: CALmod.o PPLmod.o VC.o
	$(CC) $^ -o $@ $(LAPACK)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM) $(MOD)