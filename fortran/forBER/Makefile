# Makefile for colMat.f90 noiseExist

CC = gfortran
PROGRAM = modulePPL
CFLAGS = -O3
OBJS = modulePPL.o

all: $(PROGRAM)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)