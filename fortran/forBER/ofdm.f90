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
    integer M_weight(M_tapN)
    integer M_tap(M_tapN)
    integer Mseq(Mseq_L)
    complex(kind(0d0)) MS(Nsybl,Mseq_L+1)
    complex(kind(0d0)) Tx(Mseq_L+1,GI_L+Nsybl)
    complex(kind(0d0)) Tx2((Mseq_L+1)*(GI_L+Nsybl))
    complex(kind(0d0)) Rx(Mseq_L+1+GI_L)

    !-- initialization --------------------------
    S=(0.0d0,0.0d0)
    M_tap=(/1,0,0,0/) !set M-sequence tap
    M_weight=(/1,0,0,1/) !set M-sequence weight
    Mseq=0
    MS=(0.0d0,0.0d0)

    !-- implementation --------------------------
    !M-sequence generate-------------------------
    call MseqGen(M_tap,M_weight,M_tapN,Mseq)
    !--------------------------------------------

    !set channel --------------------------------


    !--------------------------------------------

    !set data -----------------------------------
    call SetQPSKData(S,Nsybl)
    !--------------------------------------------

    !psk modulation -----------------------------
    call QPSK(S,Nsybl)
    !--------------------------------------------

    !add channel estimate sequence --------------
    call ConvStoP(Mseq,S,MS,Mseq_L,Nsybl)
    !--------------------------------------------

    !multiplexing(IFFT and GI addition) ---------
    call OFDM_Tx(MS,Tx,Nsybl,Mseq_L,GI_L)
    !--------------------------------------------

    !parallel to serial -------------------------
    call ConvPtoS()

    !--------------------------------------------

    !pass channel -------------------------------


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