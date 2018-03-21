# Makefile for myppl.f90

CC = gfortran
PROGRAM = nNoiseSeqReal nNoiseSeqMat
CFLAGS = -Wall -O3
OBJS = nNoiseSeqReal.o nNoiseSeqMat.o

all: $(PROGRAM)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< 

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: real
real:
	./nNoiseSeqReal < input.txt

.PHONY: mat
mat:
	./nNoiseSeqMat < input.txt
