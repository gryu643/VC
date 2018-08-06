# Makefile for colMat.f90 noiseExist

CC = gfortran
PROGRAM = CALmod PPLmod usePPLmod
CFLAGS = -O3
OBJS = CALmod.o PPLmod.o usePPLmod.o
MOD = pplmod.mod calmod.mod
LAPACK = -llapack

.PHONY: usePPLmod
usePPLmod: CALmod.o PPLmod.o usePPLmod.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: vc
vc: CALmod.o PPLmod.o VC.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: normal
normal: CALmod.o verify_normal.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: pdf
pdf: CALmod.o PPLmod.o eigenvaluePDF.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: pdf_ideal
pdf_ideal: eigenvaluePDF_ideal.o
	$(CC) $^ -o $@ $(LAPACK)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM) $(MOD)