# Makefile for colMat.f90 noiseExist

CC = gfortran
PROGRAM = CALmod PPLmod usePPLmod
CFLAGS = -O3
OBJS = CALmod.o PPLmod.o usePPLmod.o
MOD = calmod.mod calmod.mod

usePPLmod: CALmod.o PPLmod.o usePPLmod.o
	$(CC) $^ -o $@

VC: CALmod.o PPLmod.o VC.o
	$(CC) $^ -o $@

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)