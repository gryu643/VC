# Makefile for colMat.f90 noiseExist

CC = gfortran
PROGRAM = PPLmod usePPLmod
CFLAGS = -O3
OBJS = PPLmod.o usePPLmod.o

all: $(PROGRAM)
usePPLmod: PPLmod.o usePPLmod.o
	$(CC) $^ -o $@

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)