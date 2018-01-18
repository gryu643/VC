# Makefile for colMat.f90 noiseExist

CC = gfortran
PROGRAM = colMat
CFLAGS = -Wall -O3
OBJS = colMat.o

all: $(PROGRAM)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $<

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: mat
mat:
	./colMat < input.txt
