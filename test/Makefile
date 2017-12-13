# Makefile for myppl.f90

CC = gfortran
PROGRAM = myppl
CFLAGS = -Wall -O3
OBJS = myppl.o

$(PROGRAM): myppl.o
	$(CC) -o $@ $(CFLAGS) $<

myppl.o: myppl.f90
	$(CC) -c $<

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)