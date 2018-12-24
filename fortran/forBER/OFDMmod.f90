module OFDMmod
    use CALmod
    implicit none
contains
    subroutine SetQPSKData(S,Nsybl)
        !--------------------------------------------------------------------------!
        !set 1 or 0 to array S(QPSK case)
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) S(Nsybl)

        !-- decralation
        integer i

        !-- initialization
        S=(0.0d0,0.0d0)

        !-- implementation
        do i=1, Nsybl
            S(i) = cmplx(nint(rand()),nint(rand()),kind(0d0))
        end do
    end subroutine

    subroutine SetBPSKData(S,Nsybl)
        !--------------------------------------------------------------------------!
        !set 1 or 0 to array S(BPSK case)
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) S(Nsybl)

        !-- decralation
        integer i

        !-- initialization
        S=(0.0d0,0.0d0)

        !-- implementation
        do i=1, Nsybl
            S(i) = cmplx(nint(rand()),0.0d0,kind(0d0))
        end do
    end subroutine

    subroutine QPSK(S,Smod,Nsybl)
        !--------------------------------------------------------------------------!
        !do QPSK modulation to S set by "SetQPSKData" subroutine
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) S(Nsybl)
        complex(kind(0d0)) Smod(Nsybl)

        !-- decralation
        integer i

        !-- initialization
        Smod=(0.0d0,0.0d0)

        !-- implementation
        do i=1, Nsybl
            Smod(i) = cmplx(real(S(i))*2.0d0-1.0d0,aimag(S(i))*2.0d0-1.0d0,kind(0d0))
        end do
    end subroutine

    subroutine BPSK(S,Smod,Nsybl)
        !--------------------------------------------------------------------------!
        !do BPSK modulation to S set by "SetBPSKData" subroutine
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) S(Nsybl)
        complex(kind(0d0)) Smod(Nsybl)

        !-- decralation
        integer i

        !-- initialization
        Smod=(0.0d0,0.0d0)

        !-- implementation
        do i=1, Nsybl
            Smod(i) = cmplx(real(S(i))*2.0d0-1.0d0,0.0d0,kind(0d0))
        end do
    end subroutine

    subroutine MseqGen(M_tap,M_weight,M_tapN,Mseq)
        !--------------------------------------------------------------------------!
        !set M-sequence to array "Mseq"
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer M_tapN
        integer M_tap(M_tapN)
        integer M_weight(M_tapN)
        integer Mseq(2**M_tapN-1)
        integer M_cal
        integer TMP

        !-- decralation
        integer i,j

        !-- initialization
        TMP=0
        M_cal=0

        !-- implementation
        do i=1, 2**M_tapN-1
            !calculate xor
            M_cal = M_tap(M_tapN)*M_weight(M_tapN)
            do j=M_tapN-1, 1, -1
                M_cal = abs(M_cal-M_tap(j)*M_weight(j))
            end do

            !update M_tap
            do j=1, M_tapN
                TMP = M_tap(j)
                M_tap(j) = M_cal
                M_cal = TMP
            end do

            !output M-sequence
            Mseq(i) =  TMP
            if(Mseq(i)==0) Mseq(i)=-1
        end do
    end subroutine MseqGen

    subroutine AddMseq(Mseq,S,MS,Mseq_L,Nsybl)
        !--------------------------------------------------------------------------!
        !transform serial to parallel and add M-sequence
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Mseq_L,Nsybl
        integer Mseq(Mseq_L)
        complex(kind(0d0)) S(Nsybl)
        complex(kind(0d0)) MS(Nsybl,Mseq_L+1)

        !-- decralation
        integer i,j

        !-- initialization

        !-- implementation
        do i=1, Nsybl
            do j=1, Mseq_L
                MS(i,j) = dble(Mseq(j))
            end do
            MS(i,Mseq_L+1) = S(i)
        end do
    end subroutine AddMseq

    subroutine RemoveMseq(FFTout,MseqRx,SmodRx,Nsybl,Mseq_L)
        !--------------------------------------------------------------------------!
        !separate M-sequence and Symbol modulated
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Mseq_L,Nsybl
        complex(kind(0d0)) FFTout(Nsybl,Mseq_L+1)
        complex(kind(0d0)) MseqRx(Nsybl,Mseq_L)
        complex(kind(0d0)) SmodRx(Nsybl)

        !-- decralation
        integer i,j

        !-- initialization
        MseqRx=(0.0d0,0.0d0)
        SmodRx=(0.0d0,0.0d0)

        !-- implementation
        !store channel estimate sequence
        do j=1, Nsybl
            do i=1, Mseq_L
                MseqRx(j,i) = FFTout(j,i)
            end do
        end do

        !store Smod
        do i=1, Nsybl
            SmodRx(i) = FFTout(i,Mseq_L+1)
        end do
    end subroutine RemoveMseq

    subroutine OFDM_Tx(MS,Tx,Nsybl,Mseq_L,GI_L)
        !--------------------------------------------------------------------------!
        !inverse fast fourier transform
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,Mseq_L,GI_L
        complex(kind(0d0)) MS(Nsybl,Mseq_L+1)
        complex(kind(0d0)) Tx(Mseq_L+1,GI_L+Nsybl)

        !-- decralation
        integer i,j,k
        complex(kind(0d0)) MS_IFT(Mseq_L+1,Nsybl) !column=sample,row block

        !-- initialization
        MS_IFT=(0.0d0,0.0d0)
        Tx=(0.0d0,0.0d0)

        !-- implementation
        !IFFT
        do j=1, Mseq_L+1
            do i=1, Nsybl
                do k=1, Nsybl
                    MS_IFT(j,i) = MS_IFT(j,i) + MS(k,j) * euler(2.0d0*pi*(dble(k)/dble(Nsybl))*dble(i-1))
                end do
            end do
        end do

        !add GI
        do j=1, Mseq_L+1
            do i=1, GI_L
                Tx(j,i) = MS_IFT(j,Nsybl-GI_L+i)
            end do
            do i=1, Nsybl
                Tx(j,GI_L+i) = MS_IFT(j,i)
            end do
        end do
    end subroutine OFDM_Tx

    subroutine OFDM_Rx(Rx2,FFTout,Nsybl,Npath,Mseq_L,GI_L)
        !--------------------------------------------------------------------------!
        !fast fourier transform and GI removal
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,Npath,Mseq_L,GI_L
        complex(kind(0d0)) Rx2(Mseq_L+1,GI_L+Nsybl)
        complex(kind(0d0)) FFTout(Nsybl,Mseq_L+1)

        !-- decralation
        integer i,j,k
        complex(kind(0d0)) Rx2rmGI(Mseq_L+1,Nsybl)

        !-- initialization
        Rx2rmGI=(0.0d0,0.0d0)
        FFTout=(0.0d0,0.0d0)

        !-- implementation
        !remove GI
        do j=1, Mseq_L+1
            do i=1, Nsybl
                Rx2rmGI(j,i) = Rx2(j,GI_L+i)
            end do
        end do

        !FFT
        do j=1, Mseq_L+1
            do i=1, Nsybl
                do k=1, Nsybl
                    FFTout(i,j) = FFTout(i,j) + Rx2rmGI(j,k) * euler(-2.0d0*pi*(dble(i)/dble(Nsybl))*dble(k)) / dble(Nsybl)
                end do
            end do
        end do
    end subroutine OFDM_Rx

    subroutine ConvPtoS(Tx,Tx2,Nsybl,Mseq_L,GI_L)
        !--------------------------------------------------------------------------!
        !Convet parallel to serial
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,Mseq_L,GI_L
        complex(kind(0d0)) Tx(Mseq_L+1,GI_L+Nsybl)
        complex(kind(0d0)) Tx2((Mseq_L+1)*(GI_L+Nsybl))

        !-- decralation
        integer i,j

        !-- initialization

        !-- implementation
        do i=1, Mseq_L+1
            do j=1, GI_L+Nsybl
                Tx2((GI_L+Nsybl)*(i-1)+j) = Tx(i,j)
            end do
        end do
    end subroutine ConvPtoS

    subroutine ConvStoP(Rx,Rx2,Nsybl,Npath,Mseq_L,GI_L)
        !--------------------------------------------------------------------------!
        !Convet serial to parallel
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,Npath,Mseq_L,GI_L
        complex(kind(0d0)) Rx((Mseq_L+1)*(GI_L+Nsybl)+Npath-1)
        complex(kind(0d0)) Rx2(Mseq_L+1,GI_L+Nsybl)

        !-- decralation
        integer i,j

        !-- initialization
        Rx2=(0.0d0,0.0d0)

        !-- implementation
        do i=1, Mseq_L+1
            do j=1, GI_L+Nsybl
                Rx2(i,j) = Rx((GI_L+Nsybl)*(i-1)+j)
            end do
        end do
    end subroutine ConvStoP

    subroutine Ch(Tx2,H_weight,Rx,Nsybl,Npath,Mseq_L,GI_L)
        !--------------------------------------------------------------------------!
        !pass fading channel
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,Npath,Mseq_L,GI_L
        complex(kind(0d0)) H_weight(Npath)
        complex(kind(0d0)) Rx((Mseq_L+1)*(GI_L+Nsybl)+Npath-1)
        complex(kind(0d0)) Tx2((Mseq_L+1)*(GI_L+Nsybl))

        !-- decralation
        integer i,j
        complex(kind(0d0)) H_tap(Npath)
        integer Tx2_L

        !-- initialization
        H_tap=(0.0d0,0.0d0)
        Rx=(0.0d0,0.0d0)
        Tx2_L=(Mseq_L+1)*(GI_L+Nsybl)

        !-- implementation
        H_tap=(0.0d0,0.0d0)
        do i=1, Tx2_L+Npath-1
            !update the tap of transversal fileter
            do j=Npath, 2, -1
                H_tap(j) = H_tap(j-1)
            end do
            if(i<=Tx2_L) H_tap(1)=Tx2(i)
            if(i>Tx2_L) H_tap(1)=cmplx(0.0d0,0.0d0,kind(0d0))

            !calculate output
            do j=1, Npath
                Rx(i) = Rx(i) + H_tap(j)*H_weight(j)
            end do
        end do
    end subroutine Ch

    subroutine MakeNoise(Noise,Nsybl,Npath,Mseq_L,GI_L,KEbN0)
        !--------------------------------------------------------------------------!
        !Convet parallel to serial
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,Npath,Mseq_L,GI_L,KEbN0
        complex(kind(0d0)) Noise((Mseq_L+1)*(GI_L+Nsybl)+Npath-1)

        !-- decralation
        integer i
        integer Rx_L

        !-- initialization
        Noise=(0.0d0,0.0d0)
        Rx_L=(Mseq_L+1)*(GI_L+Nsybl)+Npath-1

        !-- implementation
        !make noise
        do i=1, Rx_L
            Noise(i) = cmplx(normal(),normal(),kind(0d0))
        end do

        !set noise by KEbN0
        Noise = Noise / sqrt(2.0d0) * sqrt(1.0d0/(10.0d0**(KEbN0/10.0d0))/2.0d0)
    end subroutine MakeNoise

    subroutine ChEstimate(MseqRx,Mseq,Chseq,Nsybl,Mseq_L)
        !--------------------------------------------------------------------------!
        !estimate channel information
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,Mseq_L
        complex(kind(0d0)) MseqRx(Nsybl,Mseq_L)
        integer Mseq(Mseq_L)
        complex(kind(0d0)) Chseq(Nsybl)

        !-- decralation
        integer i,j

        !-- initialization
        Chseq=(0.0d0,0.0d0)

        !-- implementation
        do j=1, Nsybl
            do i=1, Mseq_L
                Chseq(j) = Chseq(j) + MseqRx(j,i)*dble(Mseq(i))/dble(Mseq_L)
            end do
        end do
    end subroutine ChEstimate

    subroutine PhaseComp(Chseq,SmodRx,Nsybl)
        !--------------------------------------------------------------------------!
        !compensate the deviation of phase
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) Chseq(Nsybl)
        complex(kind(0d0)) SmodRx(Nsybl)

        !-- decralation
        integer i

        !-- initialization

        !-- implementation
        do i=1, Nsybl
            SmodRx(i) = conjg(Chseq(i))*SmodRx(i)
        end do
    end subroutine PhaseComp

    subroutine QPSKDem(SmodRx,SRx,Nsybl)
        !--------------------------------------------------------------------------!
        !QPSK demodulation
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) SmodRx(Nsybl)
        complex(kind(0d0)) SRx(Nsybl)

        !-- decralation
        integer i
        double precision Re,Im

        !-- initialization
        SRx=(0.0d0,0.0d0)

        !-- implementation
        do i=1, Nsybl
            Re=0.0d0
            Im=0.0d0
            !real part judge
            if(real(SmodRx(i))>0.0d0) then
                Re = 1.0d0
            elseif(real(SmodRx(i))>0.0d0) then
                Re = 0.0d0
            endif
            
            !imaginaly part judge
            if(aimag(SmodRx(i))>0.0d0) then
                Im = 1.0d0
            elseif(aimag(SmodRx(i))>0.0d0) then
                Im = 0.0d0
            endif

            !store Re and Im
            SRx(i) = cmplx(Re,Im,kind(0d0))
        end do
    end subroutine QPSKDem

    function cal_ber(S,SRx,Nsybl)
        !--------------------------------------------------------------------------!
        !return ber
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) S(Nsybl)
        complex(kind(0d0)) SRx(Nsybl)

        !-- decralation
        integer i,Correct,Incorrect
        double precision BER
        double precision cal_ber

        !-- initialization
        BER=0.0d0
        cal_ber=0.0d0
        Correct=0
        Incorrect=0

        !-- implementation
        do i=1, Nsybl
            !real part
            if(real(S(i))==real(SRx(i))) then
                Correct = Correct + 1
            else
                Incorrect = Incorrect + 1
            endif
            !imaginaly part
            if(aimag(S(i))==aimag(SRx(i))) then
                Correct = Correct + 1
            else
                Incorrect = Incorrect + 1
            endif
        end do

        !calculate ber
        BER = dble(Incorrect) / dble(Correct+Incorrect)
        cal_ber = BER

    end function cal_ber

    subroutine PowerAdjust(Tx2,Nsybl,Mseq_L,GI_L)
        !--------------------------------------------------------------------------!
        !adjust power(power=1)
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,Mseq_L,GI_L,i
        complex(kind(0d0)) Tx2((Mseq_L+1)*(GI_L+Nsybl))

        !-- decralation
        double precision Pow

        !-- initialization
        Pow=0.0d0

        !-- implementation
        do i=1, (Mseq_L+1)*(GI_L+Nsybl)
            Pow = Pow + abs(Tx2(i))**2 / dble((Mseq_L+1)*(GI_L+Nsybl))
        end do
        Tx2 = Tx2 / sqrt(Pow)
    end subroutine PowerAdjust

    subroutine CalPower(Rx,Noise,P,Nsybl,Npath,Mseq_L,GI_L)
        !--------------------------------------------------------------------------!
        !calculate EbN0
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,Npath,Mseq_L,GI_L,i
        complex(kind(0d0)) Rx((Mseq_L+1)*(GI_L+Nsybl)+Npath-1)
        complex(kind(0d0)) Noise((Mseq_L+1)*(GI_L+Nsybl)+Npath-1)
        double precision P(2)

        !-- decralation
        double precision Ps
        double precision Pn

        !-- initialization
        Ps=0.0d0
        Pn=0.0d0

        !-- implementation
        do i=1, (Mseq_L+1)*(GI_L+Nsybl)+Npath-1
            Ps = Ps + abs(Rx(i))**2
            Pn = Pn + abs(Noise(i))**2
        end do
        P(1) = P(1) + Ps
        P(2) = P(2) + Pn
    end subroutine CalPower
end module OFDMmod