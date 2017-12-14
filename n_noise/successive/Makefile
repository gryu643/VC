# Makefile for ppl_suc.f90

CC = gfortran
PROGRAM = ppl_suc
CFLAGS = -Wall -O3
OBJS = ppl_suc.o

$(PROGRAM): ppl_suc.o
	$(CC) -o $@ $(CFLAGS) $<

ppl_suc.o: ppl_suc.f90
	$(CC) -c $<

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: exe
exe:
	./$(PROGRAM) < input.txt
