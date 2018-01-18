# Makefile for myppl.f90

CC = gfortran
PROGRAM = seqReal seqMat
CFLAGS = -Wall -O3
OBJS = seqReal.o seqMat.o

all: $(PROGRAM)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< 

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: real
real:
	./seqReal < input.txt

.PHONY: mat
mat:
	./seqMat < input.txt
