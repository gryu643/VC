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

    subroutine QPSK(S,Nsybl)
        !--------------------------------------------------------------------------!
        !do QPSK modulation to S set by "SetQPSKData" subroutine
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) S(Nsybl)

        !-- decralation
        integer i

        !-- initialization

        !-- implementation
        do i=1, Nsybl
            S(i) = cmplx(real(S(i))*2.0d0-1.0d0,aimag(S(i))*2.0d0-1.0d0,kind(0d0))
        end do
    end subroutine

    subroutine BPSK(S,Nsybl)
        !--------------------------------------------------------------------------!
        !do BPSK modulation to S set by "SetBPSKData" subroutine
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) S(Nsybl)

        !-- decralation
        integer i

        !-- initialization

        !-- implementation
        do i=1, Nsybl
            S(i) = cmplx(real(S(i))*2.0d0-1.0d0,0.0d0,kind(0d0))
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

    subroutine ConvStoP(Mseq,S,MS,Mseq_L,Nsybl)
        !--------------------------------------------------------------------------!
        !transform serial to parallel
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
    end subroutine ConvStoP

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
                    MS_IFT(j,i) = MS_IFT(j,i) + MS(k,j) * euler(2.0d0*pi*(dble(k)/dble(Nsybl))*i)
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
end module OFDMmod