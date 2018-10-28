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
    integer,parameter :: Nloop=10000
    double precision,parameter :: BERStandard=1.0d-2

    integer i,j
    double precision Ampd(Npath,1)
    double precision Psig
    double precision Pwgn
    integer Collect
    integer False
    integer loop
    integer KEbN0
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
    complex(kind(0d0)) Xppl(Nsybl,Nsybl)
    complex(kind(0d0)) HH(Nsybl,Nsybl+Npath-1)
    complex(kind(0d0)) HHH(Nsybl,Nsybl)
    complex(kind(0d0)) S(1,Nsybl)
    integer TdatI(1,Nsybl)
    integer TdatQ(1,Nsybl)
    complex(kind(0d0)) V(Nsybl,Nsybl)
    complex(kind(0d0)) SU(Nsybl,Nsybl)
    complex(kind(0d0)) X(Nsybl,Nsybl)
    double precision Pow
    complex(kind(0d0)) Y(Nsybl+Npath-1,1)
    complex(kind(0d0)) Noise(Nsybl+Npath-1,1)
    complex(kind(0d0)) Yn(Nsybl,1)
    complex(kind(0d0)) Y2(Nsybl,1)
    integer RdatI(1,Nsybl)
    integer RdatQ(1,Nsybl)
    complex(kind(0d0)) A(Nsybl,1)
    complex(kind(0d0)) Y2H(1,Nsybl)
    complex(kind(0d0)) R(1,1)
    complex(kind(0d0)) R2
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
    Psig=0.0d0
    Pwgn=0.0d0
    Collect=0
    False=0
    loop=0
    KEbN0=0
    Cpath=(0.0d0,0.0d0)
    H=(0.0d0,0.0d0)
    HE=(0.0d0,0.0d0)
    Xppl=(0.0d0,0.0d0)
    HH=(0.0d0,0.0d0)
    HHH=(0.0d0,0.0d0)
    S=0.0d0
    TdatI=0
    TdatQ=0
    V=(0.0d0,0.0d0)
    SU=(0.0d0,0.0d0)
    X=(0.0d0,0.0d0)
    Pow=0.0d0
    Y=(0.0d0,0.0d0)
    Noise=(0.0d0,0.0d0)
    Yn=(0.0d0,0.0d0)
    Y2=(0.0d0,0.0d0)
    RdatI=0
    RdatQ=0
    A=(0.0d0,0.0d0)
    Y2H=(0.0d0,0.0d0)
    R=(0.0d0,0.0d0)
    R2=(0.0d0,0.0d0)
    EbN0=0.0d0
    BER=0.0d0
    Eig=0.0d0
    InstantBER=0.0d0
    LambdaEbN0=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open (1, file='VCQ(s32p4).csv', status='replace')
    open (2, file='VCtrial.csv', status='replace')
    open (3, file='VCQ_UseChNum.csv', status='replace')

    !implimentation part
    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do

    do KEbN0=SEbN0, EEbN0, Step !Eb/N0 loop
        Psig = 0.0d0
        Pwgn = 0.0d0
        Collect = 0
        False = 0
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

            !set information symbol
            do i=1, Nsybl
                S(1,i) = cmplx(nint(rand())*2.0d0-1.0d0,nint(rand())*2.0d0-1.0d0,kind(0d0)) !-1or1
            end do
            do i=1, Nsybl
                TdatI(1,i) = (real(S(1,i))+1)/2 !0or1
                TdatQ(1,i) = (aimag(S(1,i))+1)/2
            end do
            do i=1, Nsybl
                do j=1, Nsybl
                    SU(j,i) = S(1,i) * V(j,i)
                end do
            end do

            !transmit vector
            X = (0.0d0,0.0d0)
            do i=1,Nsybl
                do j=1,Nsybl
                    X(j,1) = X(j,1) + SU(j,i)
                end do
            end do

            !power
            Pow = 0.0d0
            do i=1, Nsybl
                Pow = Pow + (abs(X(i,1))**2.0d0)/Nsybl
            end do
            X = X / sqrt(Pow)
            call CMultiply(H,X,Y,Nsybl+Npath-1,Nsybl,Nsybl,1) !pass H

            do i=1, Nsybl+Npath-1
                Noise(i,1) = cmplx(normal(),normal(),kind(0d0))
            end do
            !正規乱数の分散＝１＝雑音電力なので、正規乱数に雑音電力をかける（√２で割っているのはIとQの両方合わせて雑音電力とするため）
            Noise = Noise / sqrt(2.0d0) * sqrt(1.0d0/(10.0d0**(KEbN0/10.0d0))/2.0d0)

            do i=1, Nsybl+Npath-1
                Psig = Psig + (abs(Y(i,1))**2.0d0)
                Pwgn = Pwgn + abs(Noise(i,1))**2.0d0
            end do

            !add noise
            call CAdd(Y,Noise,Y,Nsybl+Npath-1,1,Nsybl+Npath-1,1)

            !Matched Filter
            call CMultiply(HH,Y,Yn,Nsybl,Nsybl+Npath-1,Nsybl+Npath-1,1)
        
            do i=1, Nsybl
                Y2(i,1) = Yn(i,1)
            end do
            RdatI = 0
            RdatQ = 0
            !Demodulation
            do i=1, Nsybl
                do j=1, Nsybl
                    A(j,1) = V(j,i)
                end do
                !inner product
                call CAdjoint(Y2,Y2H,Nsybl,1)
                call CMultiply(Y2H,A,R,1,Nsybl,Nsybl,1)
                R(1,1) = conjg(R(1,1))
                R2 = R(1,1)

                if(real(R2)>0.0) then
                    RdatI(1,i) = 1
                elseif(real(R2)<0.0) then
                    RdatI(1,i) = 0
                endif
                if(aimag(R2)>0.0) then
                    RdatQ(1,i) = 1
                elseif(aimag(R2)<0.0) then
                    RdatQ(1,i) = 0
                endif
            end do

            !Bit Error Rate
            do i=1, Nsybl
                if(RdatI(1,i)==TdatI(1,i)) then
                    Collect = Collect + 1
                elseif(RdatI(1,i)/=TdatI(1,i)) then
                    False = False + 1
                endif
                if(RdatQ(1,i)==TdatQ(1,i)) then
                    Collect = Collect + 1
                elseif(RdatQ(1,i)/=TdatQ(1,i)) then
                    False = False + 1
                endif
            end do
            write(2,*) loop, ',', dble(False)/(dble(Collect)+dble(False))
        end do
        
        EbN0 = 10.0d0*dlog10(Psig/Pwgn/2.0d0) !QPSK rate =2
        BER = dble(False) / (dble(Collect) + dble(False))
        write(1,*) EbN0, ',', BER, ',', RTNum
        print *, 'EbN0=', EbN0
        print *, 'BER=', BER
        print *, 'AvRTNum=', AvRTNum
        write(3,*) KEbN0, ',', AvUseChNum
    end do

    !file close
    close(1)
    close(2)
    close(3)

    !time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program VCQ