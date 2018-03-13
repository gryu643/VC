# Makefile for colMat.f90 noiseExist

CC = gfortran
PROGRAM = eNoiseColMat
CFLAGS = -Wall -O3
OBJS = eNoiseColMat.o

all: $(PROGRAM)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $<

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: run
run:
	./eNoiseColMat < input.txt
