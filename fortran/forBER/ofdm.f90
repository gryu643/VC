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

    integer i,j
    complex(kind(0d0)) S(Nsybl)
    complex(kind(0d0)) Smod(Nsybl)
    complex(kind(0d0)) Sp(Nsybl,1)
    integer M_weight(M_tapN)
    integer M_tap(M_tapN)
    integer Mseq(Mseq_L)
    complex(kind(0d0)) MS(Nsybl,Mseq_L+1)
    complex(kind(0d0)) Tx(Mseq_L+1,GI_L+Nsybl)
    complex(kind(0d0)) Tx2((Mseq_L+1)*(GI_L+Nsybl))
    complex(kind(0d0)) Rx((Mseq_L+1)*(GI_L+Nsybl)+Npath-1)
    double precision Ampd(Npath)
    complex(kind(0d0)) H_weight(Npath)
    integer KEbN0
    complex(kind(0d0)) Noise((Mseq_L+1)*(GI_L+Nsybl)+Npath-1)

    !-- initialization --------------------------
    S=(0.0d0,0.0d0)
    Smod=(0.0d0,0.0d0)
    M_tap=(/1,0,0,0/) !set M-sequence tap
    M_weight=(/1,0,0,1/) !set M-sequence weight
    Mseq=0
    MS=(0.0d0,0.0d0)
    Sp=(0.0d0,0.0d0)
    Ampd=0.0d0
    H_weight=(0.0d0,0.0d0)
    Noise=(0.0d0,0.0d0)

    !-- implementation --------------------------
    !M-sequence generate-------------------------
    call MseqGen(M_tap,M_weight,M_tapN,Mseq)
    !--------------------------------------------

    !channel gain parameter ---------------------
    do i=1, Npath
        Ampd(i) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do
    !--------------------------------------------

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

    !pass channel -------------------------------
    call Ch(Tx2,H_weight,Rx,Nsybl,Npath,Mseq_L,GI_L)
    !--------------------------------------------

    !add noise ----------------------------------
    call MakeNoise(Noise,Nsybl,Npath,Mseq_L,GI_L,KEbN0)
    Rx = Rx + Noise
    !--------------------------------------------

    !serial to parallel -------------------------

    !--------------------------------------------

    !remove GI ----------------------------------


    !--------------------------------------------

    !FFT ----------------------------------------


    !--------------------------------------------

    !remove channel estimate sequence -----------


    !--------------------------------------------

    !channel estimate ---------------------------


    !--------------------------------------------

    !phase compensation -------------------------


    !--------------------------------------------

    !parallel to serial -------------------------


    !--------------------------------------------

    !PSK demodulation ---------------------------


    !--------------------------------------------
end program ofdm