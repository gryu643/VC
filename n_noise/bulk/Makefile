# Makefile for myppl.f90

CC = gfortran
PROGRAM = ppl_bulk ppl_abs_bulk
CFLAGS = -Wall -O3
OBJS = ppl_bulk.o ppl_abs_bulk.o

ppl_bulk: ppl_bulk.o
	$(CC) -o $@ $<

ppl_abs_bulk: ppl_abs_bulk.o
	$(CC) -o $@ $<

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $<

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: exe
exe:
	./$(PROGRAM) < input.txt
