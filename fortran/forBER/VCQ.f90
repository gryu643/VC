program VCQ
    use PPLVCQmod
    use CALmod
    implicit none

    !run time declaration
    integer t1, t2, t_rate, t_max, diff

    !declaration
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=4
    integer,parameter :: SEbN0=-10
    integer,parameter :: EEbN0=40
    integer,parameter :: Step=5
    integer,parameter :: Nloop=10
    double precision,parameter :: BERStandard=1.0d-2

    integer i,j
    double precision Ampd(Npath,1)
    integer loop
    integer KEbN0
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
    complex(kind(0d0)) Xppl(Nsybl,Nsybl)
    complex(kind(0d0)) HH(Nsybl,Nsybl+Npath-1)
    complex(kind(0d0)) HHH(Nsybl,Nsybl)
    complex(kind(0d0)) V(Nsybl,Nsybl)
    double precision EbN0
    double precision BER
    double precision Eig(1,Nsybl)
    double precision BER_TMP
    double precision AvUseChNum
    double precision InstantBER
    double precision LambdaEbN0
    double precision ConvStandard
    integer RTNum
    double precision AvRTNum

    !initialize
    Ampd=0.0d0
    loop=0
    KEbN0=0
    Cpath=(0.0d0,0.0d0)
    H=(0.0d0,0.0d0)
    HE=(0.0d0,0.0d0)
    Xppl=(0.0d0,0.0d0)
    HH=(0.0d0,0.0d0)
    HHH=(0.0d0,0.0d0)
    V=(0.0d0,0.0d0)
    EbN0=0.0d0
    BER=0.0d0
    Eig=0.0d0
    InstantBER=0.0d0
    LambdaEbN0=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open (3, file='VCQ_UseChNum.csv', status='replace')

    !implimentation part
    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do

    do KEbN0=SEbN0, EEbN0, Step !Eb/N0 loop
        AvUseChNum=0.0d0
        AvRTNum=0.0d0
        !setup Convergence standard
        if(KEbN0<0) then
            ConvStandard = 1.0d0
        else
            ConvStandard = 1.0d0*dexp(-0.23*dble(KEbN0))
        endif

        do loop=1, Nloop !Monte calro loop
            do i=1, Npath
                Cpath(i,1) = cmplx(normal(),normal(),kind(0d0))/sqrt(2.0d0)*Ampd(i,1)
            end do

            !set H
            do i=1, Nsybl
                do j=1, Npath
                    H(i+j-1,i) = Cpath(j,1)
                end do
            end do

            !set HE
            do i=1, Nsybl+Npath-1
                do j=1, Npath
                    HE(i+j-1,i) = Cpath(j,1)
                end do
            end do

            !set Xppl
            do i=1, Nsybl
                do j=1, Nsybl
                    Xppl(i,j) = cmplx(1.0d0, 0.0d0, kind(0d0))
                end do
            end do

            !set HH
            call CAdjoint(H,HH,Nsybl+Npath-1,Nsybl)

            !set HHH
            call CMultiply(HH,H,HHH,Nsybl,Nsybl+Npath-1,Nsybl+Npath-1,Nsybl)

            !eigenvalue decomposition
            call PPLVCQ(H,HE,Xppl,Eig,Nsybl,Npath,RTNum,ConvStandard)
            call CSubstitute(V,Xppl,Nsybl,Nsybl)

            AvRTNum = AvRTNum + dble(RTNum)/dble(Nloop)

            !judge quality
            do i=1, Nsybl
                BER_TMP=0.0d0
                do j=1, i
                    LambdaEbN0 = Eig(1,j)*10.0d0**(dble(KEbN0)/10.0d0)
                    InstantBER = 1.0d0/2.0d0*erfc(sqrt(LambdaEbN0))
                    BER_TMP = BER_TMP + InstantBER/dble(i)
                end do

                if(BER_TMP>BERStandard) then
                    AvUseChNum = AvUseChNum + dble(i-1)/dble(Nloop)
                    exit
                endif
                
                if(i==Nsybl) then
                    AvUseChNum = AvUseChNum + dble(Nsybl)/dble(Nloop)
                endif
            end do
        end do

        print *, 'KEbN0=', KEbN0
        print *, 'AvRTNum=', AvRTNum
        write(3,*) KEbN0, ',', AvUseChNum
    end do

    !file close
    close(3)

    !time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program VCQ