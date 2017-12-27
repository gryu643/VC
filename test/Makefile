# Makefile for myppl.f90

CC = gfortran
PROGRAM = myppl
CFLAGS = -O3 -ffree-form
OBJS = myppl.o

$(PROGRAM): myppl.o
	$(CC) -o $@ $(CFLAGS) $<

myppl.o: myppl.f90
	$(CC) -c $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM)

.PHONY: exe
exe:
	./$(PROGRAM) < input.txt
