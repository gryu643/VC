# Makefile for myppl.f90

CC = gfortran
PROGRAM = colReal colMat
CFLAGS = -Wall -O3
OBJS = colReal.o colMat.o

all: $(PROGRAM)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $<

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: real
real:
	./colReal < input.txt

.PHONY: mat
mat:
	./colMat < input.txt
