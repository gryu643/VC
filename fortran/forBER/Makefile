# Makefile for colMat.f90 noiseExist

CC = gfortran
PROGRAM = CALmod PPLmod usePPLmod vc vc_nonfading normal pdf pdf_ideal qpsk
CFLAGS = -O3
OBJS = *.o
MOD = *.mod
LAPACK = -llapack -lblas

.PHONY: usePPLmod
usePPLmod: CALmod.o PPLmod.o usePPLmod.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: vc
vc: CALmod.o PPLmod.o VC.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: vc_nonfading
vc_nonfading: CALmod.o PPLmod.o VC_nonfading.o
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

.PHONY: qpsk
qpsk: CALmod.o qpsk.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: qpsk_ideal
qpsk_ideal: qpsk_ideal.o
	$(CC) $^ -o $@ $(LAPACK)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM) $(MOD)