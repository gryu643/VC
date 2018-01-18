# Makefile for myppl.f90

CC = gfortran
PROGRAM = ppl_suc ppl_abs_suc
CFLAGS = -Wall -O3
OBJS = ppl_suc.o ppl_abs_suc.o

all: $(PROGRAM)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< 

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: real
real:
	./ppl_suc < input.txt

.PHONY: mat
mat:
	./ppl_abs_suc < input.txt