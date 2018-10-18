# Makefile for colMat.f90 noiseExist

CC = gfortran
PROGRAM = CALmod PPLmod usePPLmod vc vc_nonfading normal pdf pdf_ideal qpsk qpsk_ideal_fading vc_ideal vc_ideal_formula
CFLAGS = -O3 -mtune=native -march=native -mfpmath=both
OBJS = *.o
MOD = *.mod
LAPACK = -llapack -lblas

.PHONY: usePPLmod
usePPLmod: CALmod.o PPLCUTOFFmod.o usePPLmod.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: avoth_ber
avoth_ber: CALmod.o PPLCUTOFFmod.o AvOth_BER.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: vc
vc: CALmod.o PPLmod.o VC.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: cutoff
cutoff: CALmod.o PPLCutoffmod.o VC_cutoff.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: vc_ideal
vc_ideal: VC_ideal.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: vc_ideal2
vc_ideal2: CALmod.o PPLmod.o VC_ideal_v2.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: vc_ideal_formula
vc_ideal_formula: VC_ideal_formula.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: each_ch
each_ch: CALmod.o PPLmod.o VC_BEReachCH.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: fixedH
fixedH: CALmod.o PPLmod.o VC_fixedH.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: Pcon
Pcon: CALmod.o PPLmod.o PCONmod2.o PCONmod.o VC_Pcontrol.o
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

.PHONY: qpsk_fading
qpsk_fading: CALmod.o qpsk_fading.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: qpsk_ideal_fading
qpsk_ideal_fading: CALmod.o qpsk_ideal_fading.o
	$(CC) $^ -o $@ $(LAPACK)

.PHONY: MRdiv
MRdiv: CALmod.o qpsk_ideal_MRdiv.o
	$(CC) $^ -o $@ $(LAPACK)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM) $(MOD)