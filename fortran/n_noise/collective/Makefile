# Makefile for myppl.f90

CC = gfortran
PROGRAM = nNoiseColReal nNoiseColMat
CFLAGS = -Wall -O3
OBJS = nNoiseColReal.o nNoiseColMat.o

all: $(PROGRAM)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $<

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: real
real:
	./nNoiseColReal < input.txt

.PHONY: mat
mat:
	./nNoiseColMat < input.txt
