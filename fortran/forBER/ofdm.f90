program ofdm
    use CALmod
    use OFDMmod
    implicit none

    !-- declaration -----------------------------
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=2
    integer,parameter :: M_tapN=4
    integer,parameter :: GI_L=8
    integer,parameter :: Mseq_L=2**M_tapN-1
    integer,parameter :: Nloop=1
    integer,parameter :: SEbN0=20
    integer,parameter :: EEbN0=20
    integer,parameter :: Step=10

    integer i
    complex(kind(0d0)) S(Nsybl)
    complex(kind(0d0)) SRx(Nsybl)
    complex(kind(0d0)) Smod(Nsybl)
    complex(kind(0d0)) SmodRx(Nsybl)
    complex(kind(0d0)) Sp(Nsybl,1)
    integer M_weight(M_tapN)
    integer M_tap(M_tapN)
    integer Mseq(Mseq_L)
    complex(kind(0d0)) MseqRx(Nsybl,Mseq_L)
    complex(kind(0d0)) MS(Nsybl,Mseq_L+1)
    complex(kind(0d0)) Tx(Mseq_L+1,GI_L+Nsybl) !Tx(parallel)
    complex(kind(0d0)) Tx2((Mseq_L+1)*(GI_L+Nsybl)) !Tx(serial)
    complex(kind(0d0)) Rx((Mseq_L+1)*(GI_L+Nsybl)+Npath-1) !Rx(serial)
    complex(kind(0d0)) Rx2(Mseq_L+1,GI_L+Nsybl) !Rx(parallel)
    double precision Ampd(Npath)
    complex(kind(0d0)) H_weight(Npath)
    integer KEbN0
    integer Kloop
    complex(kind(0d0)) Noise((Mseq_L+1)*(GI_L+Nsybl)+Npath-1)
    complex(kind(0d0)) FFTout(Nsybl,Mseq_L+1) !FFT output
    complex(kind(0d0)) Chseq(Nsybl)
    double precision BER

    !-- initialization --------------------------
    S=(0.0d0,0.0d0)
    Smod=(0.0d0,0.0d0)
    SmodRx=(0.0d0,0.0d0)
    M_tap=(/1,0,0,0/) !set M-sequence tap
    M_weight=(/1,0,0,1/) !set M-sequence weight
    Mseq=0
    MS=(0.0d0,0.0d0)
    Sp=(0.0d0,0.0d0)
    Ampd=0.0d0
    H_weight=(0.0d0,0.0d0)
    Noise=(0.0d0,0.0d0)
    FFTout=(0.0d0,0.0d0)
    MseqRx=(0.0d0,0.0d0)
    Chseq=(0.0d0,0.0d0)
    KEbN0=-10
    BER=0.0d0

    !-- implementation --------------------------
    !file open ----------------------------------
    open (1,file='ofdm.csv',status='replace')
    !--------------------------------------------

    !M-sequence generate-------------------------
    call MseqGen(M_tap,M_weight,M_tapN,Mseq)
    !--------------------------------------------

    !channel gain parameter ---------------------
    do i=1, Npath
        Ampd(i) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do
    !--------------------------------------------

    do KEbN0=SEbN0, EEbN0, Step
        BER=0.0d0

        do Kloop=1, Nloop
            !set channel --------------------------------
            do i=1, Npath
                H_weight(i) = cmplx(normal(),normal(),kind(0d0))/sqrt(2.0d0)*Ampd(i)
            end do
            !--------------------------------------------

            !set data -----------------------------------
            call SetQPSKData(S,Nsybl)
            !--------------------------------------------

            !psk modulation -----------------------------
            call QPSK(S,Smod,Nsybl)
            !--------------------------------------------

            !add channel estimate sequence --------------
            call AddMseq(Mseq,Smod,MS,Mseq_L,Nsybl)
            !--------------------------------------------

            !multiplexing(IFFT and GI addition) ---------
            call OFDM_Tx(MS,Tx,Nsybl,Mseq_L,GI_L)
            !--------------------------------------------

            !parallel to serial -------------------------
            call ConvPtoS(Tx,Tx2,Nsybl,Mseq_L,GI_L)
            !--------------------------------------------

            !power=1 ------------------------------------
            call PowerAdjust(Tx2,Nsybl,Mseq_L,GI_L)
            !--------------------------------------------

            !pass channel -------------------------------
            call Ch(Tx2,H_weight,Rx,Nsybl,Npath,Mseq_L,GI_L)
            !--------------------------------------------

            !add noise ----------------------------------
            call MakeNoise(Noise,Nsybl,Npath,Mseq_L,GI_L,KEbN0)
            Rx = Rx + Noise
            !--------------------------------------------

            !serial to parallel -------------------------
            call ConvStoP(Rx,Rx2,Nsybl,Npath,Mseq_L,GI_L)
            !--------------------------------------------

            !demultiplexing(FFT and GI removal) ---------
            call OFDM_Rx(Rx2,FFTout,Nsybl,Npath,Mseq_L,GI_L)
            !--------------------------------------------

            !remove channel estimate sequence -----------
            call RemoveMseq(FFTout,MseqRx,SmodRx,Nsybl,Mseq_L)
            !--------------------------------------------

            !channel estimate ---------------------------
            call ChEstimate(MseqRx,Mseq,Chseq,Nsybl,Mseq_L)
            !--------------------------------------------

            !phase compensation -------------------------
            call PhaseComp(Chseq,SmodRx,Nsybl)
            !--------------------------------------------

            !PSK demodulation ---------------------------
            call QPSKDem(SmodRx,SRx,Nsybl)
            !--------------------------------------------

            !calculate BER ------------------------------
            BER = BER + cal_ber(S,SRx,Nsybl)/dble(Nloop)
            !--------------------------------------------
        end do

        !result output ------------------------------
        print *, KEbN0
        print *, BER
        !--------------------------------------------
    end do

    !file close ---------------------------------
    close(1)
    !--------------------------------------------
end program ofdm