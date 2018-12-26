module CombSubmod
    use CALmod
    implicit none
contains 
    subroutine CombOth(X,BERFlag,Ksybl,AVGOTH,AVGOTHN,ConvSize,Nsybl)
        !------------------------------------------------------------------!
        !calculate orthogonality(subroutine used in PPLCombmod.f90)
        !------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,ConvSize
        complex(kind(0d0)) X(Nsybl,Nsybl)
        logical BERFlag(ConvSize)
        double precision AVGOTH(ConvSize) !各雑音条件（Ksybl分）での直交度
        double precision AVGOTHN !全シンボル分での直交度
        integer Ksybl(ConvSize)

        !-- declaration
        integer i,j,k,m
        complex(kind(0d0)) V1(Nsybl),V2(Nsybl)
        complex(kind(0d0)) NAISEKI(1,1)
        double precision NAISEKI_TMP

        !-- initialization
        NAISEKI=0.0d0
        NAISEKI_TMP=0.0d0
        V1=(0.0d0,0.0d0)
        V2=(0.0d0,0.0d0)

        !-- implementation
        do i=2, Nsybl
            do j=1, i-1
                !固有ベクトル群を１列のベクトルに格納
                do k=1, Nsybl
                    V1(k) = X(k,i)
                end do
                !内積を取る固有ベクトルを格納
                do k=1, Nsybl
                    V2(k) = X(k,j)
                end do
                NAISEKI(1,1) = NAISEKI(1,1) + abs(dot_product(V1,V2))

                call CAbs(NAISEKI,NAISEKI_TMP,1,1)
                do m=1, ConvSize
                    if(BERFlag(m)) cycle

                    if(i==Ksybl(m)) then
                        AVGOTH(m) = NAISEKI_TMP
                    endif
                end do
            end do
        end do
        call CAbs(NAISEKI,NAISEKI_TMP,1,1)
        AVGOTHN = NAISEKI_TMP

    end subroutine CombOth

    subroutine CombOutMax(X,LAMBDA,BERFlag,Nsybl,UseChNum,ConvSize,l,V,Eig,Pt,RTNum)
        !-----------------------------------------------------------------------------!
        !output when l==MaxPocnLoop(subroutine used in PPLCombmod.f90)
        !-----------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,l,ConvSize
        complex(kind(0d0)) X(Nsybl,ConvSize)
        complex(kind(0d0)) LAMBDA(Nsybl,1)
        logical BERFlag(ConvSize)
        integer UseChNum(ConvSize)
        complex(kind(0d0)) V(ConvSize,Nsybl,Nsybl)
        double precision Eig(Nsybl,ConvSize)
        double precision Pt(ConvSize,Nsybl)
        integer RTNum(ConvSize)

        !-- declaration
        integer j,k,m

        !-- initialization

        !-- implementation
        do j=1, ConvSize
            if(BERFlag(j)) cycle

            BERFlag(j) = .True.
            UseChNum(j) = Nsybl
            RTNum(j) = l
            do m=1, Nsybl
                Eig(m,j) = real(LAMBDA(m,1))
            end do

            !最大反復数で打ち切る際は、送信電力分配しない
            do k=1, Nsybl
                Pt(j,k) = 1.0d0
            end do

            do m=1, Nsybl
                do k=1, Nsybl
                    V(j,k,m) = X(k,m)
                end do
            end do
        end do

    end subroutine CombOutMax

    subroutine CombOut(X,LAMBDA,RTNum,Eig,UseChNum,V,Nsybl,ConvSize,Addr,l)
        !-----------------------------------------------------------------------------!
        !output(normal case)(subroutine used in PPLCombmod.f90)
        !-----------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl,Addr,l,ConvSize
        complex(kind(0d0)) X(Nsybl,ConvSize)
        complex(kind(0d0)) LAMBDA(Nsybl,1)
        integer RTNum(ConvSize)
        double precision Eig(Nsybl,ConvSize)
        integer UseChNum(ConvSize)
        complex(kind(0d0)) V(ConvSize,Nsybl,Nsybl)

        !-- decralation
        integer k,m

        !-- initialization

        !-- implementation
        RTNum(Addr) = l
        do m=1, Nsybl
            Eig(m,Addr) = real(LAMBDA(m,1))
        end do

        if(UseChNum(Addr)/=0) then
            do m=1, UseChNum(Addr)
                do k=1, Nsybl
                    V(Addr,k,m) = X(k,m)
                end do
            end do
        endif
    end subroutine CombOut

    subroutine CombIdeal(lambda,EbN0,Pt_TMP,BER)
        !-----------------------------------------------------------------------------!
        !calculate Ideal BER(subroutine used in PPLCombmod.f90)
        !-----------------------------------------------------------------------------!
        implicit none

        !-- argument
        double precision lambda
        double precision EbN0
        double precision Pt_TMP
        double precision BER

        !-- decralation
        double precision LambdaEbN0
        double precision InstantBER

        !-- initialization
        LambdaEbN0=0.0d0
        InstantBER=0.0d0

        !-- implementation
        LambdaEbN0 = lambda*EbN0*Pt_TMP
        InstantBER = 1.0d0/2.0d0*erfc(sqrt(LambdaEbN0))
        BER = BER + InstantBER

    end subroutine CombIdeal
end module CombSubmod
