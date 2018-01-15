# Makefile for myppl.f90

CC = gfortran
PROGRAM = ppl_bulk ppl_abs_bulk
CFLAGS = -Wall -O3
OBJS = ppl_bulk.o ppl_abs_bulk.o

all: $(PROGRAM)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $<

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: exe
exe:
	./ppl_bulk < input.txt
