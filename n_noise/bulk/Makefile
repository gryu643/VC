# Makefile for myppl.f90

CC = gfortran
PROGRAM = ppl_bulk
CFLAGS = -Wall -O3
OBJS = ppl_bulk.o

$(PROGRAM): ppl_bulk.o
	$(CC) -o $@ $<

ppl_bulk.o: ppl_bulk.f90
	$(CC) -c $< 

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: exe
exe:
	./$(PROGRAM) < input.txt
