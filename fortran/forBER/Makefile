# Makefile for colMat.f90 noiseExist

CC = gfortran
PROGRAM = usePPL avoth_ber tsym vc vcq cutoff comb vc_ideal vc_ideal2 vc_ideal_formula \
each_ch fixedH pcon pcon2 vc_nonfading normal pdf pdf_ideal qpsk qpsk_ideal qpsk_fading \
rake est est_err qpsk_ideal_fading MRdiv ofdm
CFLAGS = -O3 -mtune=native -march=native -mfpmath=both
OBJS = *.o
MOD = *.mod
LAPACK = -llapack -lblas

.PHONY: $(PROGRAM)
usePPL: CALmod.o ChEstmod.o PPLNmod.o PPLmod.o usePPL.o
	$(CC) $^ -o $@ $(LAPACK)

avoth_ber: CALmod.o PPLAvOth_BERmod.o AvOth_BER.o
	$(CC) $^ -o $@ $(LAPACK)

tsym: CALmod.o PPLTSYMmod.o AVOth_BER_TSYM.o
	$(CC) $^ -o $@ $(LAPACK)

vc: CALmod.o PPLmod.o VC.o
	$(CC) $^ -o $@ $(LAPACK)

vcq: CALmod.o PPLVCQmod.o VCQ.o
	$(CC) $^ -o $@ $(LAPACK)

cutoff: CALmod.o PPLCutoffmod.o VC_cutoff.o
	$(CC) $^ -o $@ $(LAPACK)

comb: CALmod.o PCONmod2.o CombSubmod.o PPLCombmod.o VC_comb.o
	$(CC) $^ -o $@ $(LAPACK)

vc_ideal: VC_ideal.o
	$(CC) $^ -o $@ $(LAPACK)

vc_ideal2: CALmod.o PPLmod.o VC_ideal_v2.o
	$(CC) $^ -o $@ $(LAPACK)

vc_ideal_formula: VC_ideal_formula.o
	$(CC) $^ -o $@ $(LAPACK)

each_ch: CALmod.o PPLmod.o VC_BEReachCH.o
	$(CC) $^ -o $@ $(LAPACK)

fixedH: CALmod.o PPLmod.o VC_fixedH.o
	$(CC) $^ -o $@ $(LAPACK)

pcon: CALmod.o PPLmod.o PCONmod.o PCONmod2.o VC_Pcontrol.o
	$(CC) $^ -o $@ $(LAPACK)

pcon2: CALmod.o PPLmod.o PCONmod2.o PCONmod3.o VC_Pcontrol2.o
	$(CC) $^ -o $@ $(LAPACK)

vc_nonfading: CALmod.o PPLmod.o VC_nonfading.o
	$(CC) $^ -o $@ $(LAPACK)

normal: CALmod.o verify_normal.o
	$(CC) $^ -o $@ $(LAPACK)

pdf: CALmod.o PPLmod.o eigenvaluePDF.o
	$(CC) $^ -o $@ $(LAPACK)

pdf_ideal: eigenvaluePDF_ideal.o
	$(CC) $^ -o $@ $(LAPACK)

qpsk: CALmod.o qpsk.o
	$(CC) $^ -o $@ $(LAPACK)

qpsk_ideal: qpsk_ideal.o
	$(CC) $^ -o $@ $(LAPACK)

qpsk_fading: CALmod.o qpsk_fading.o
	$(CC) $^ -o $@ $(LAPACK)

rake: CALmod.o rake.o
	$(CC) $^ -o $@ $(LAPACK)

est: CALmod.o ChEstmod.o ch_est.o
	$(CC) $^ -o $@ $(LAPACK)

est_err: CALmod.o ChEstmod.o ch_est_err.o
	$(CC) $^ -o $@ $(LAPACK)

qpsk_ideal_fading: CALmod.o qpsk_ideal_fading.o
	$(CC) $^ -o $@ $(LAPACK)

MRdiv: CALmod.o qpsk_ideal_MRdiv.o
	$(CC) $^ -o $@ $(LAPACK)

ofdm: CALmod.o OFDMmod.o ofdm.o
	$(CC) $^ -o $@ $(LAPACK)

.SUFFIXES: .f90 .o
.f90.o:
	$(CC) -c $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJS) $(PROGRAM) $(MOD)