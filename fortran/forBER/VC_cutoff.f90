program AvOth_BER
    use PPLCUTOFFmod
    use CALmod
    implicit none

    !run time declaration
    integer t1, t2, t_rate, t_max, diff

    !declaration
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=4
    integer,parameter :: Nloop=100
    integer,parameter :: SEbN0=-10
    integer,parameter :: EEbN0=40
    integer,parameter :: Step=5
    double precision,parameter :: BERStandard=1.0d-2
    integer,parameter :: ConvSize=(int((EEbN0-SEbN0)/Step)+1)

    integer i,j,k
    integer RTNum(ConvSize), UseChNum(ConvSize)
    integer KEbN0
    double precision Ampd(Npath,1)
    double precision Psig(ConvSize)
    double precision Pwgn(ConvSize)
    integer Collect(ConvSize)
    integer False(ConvSize)
    integer loop
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
    complex(kind(0d0)) Xppl(Nsybl,Nsybl)
    complex(kind(0d0)) HH(Nsybl,Nsybl+Npath-1)
    complex(kind(0d0)) HHH(Nsybl,Nsybl)
    complex(kind(0d0)) S(1,Nsybl)
    integer TdatI(Nsybl)
    integer TdatQ(Nsybl)
    complex(kind(0d0)) V(Nsybl,Nsybl)
    complex(kind(0d0)) PV(ConvSize,Nsybl,Nsybl)
    complex(kind(0d0)) SU(Nsybl,Nsybl)
    complex(kind(0d0)) X(Nsybl,Nsybl)
    double precision Pow
    complex(kind(0d0)) Y(Nsybl+Npath-1,1)
    complex(kind(0d0)) Noise(Nsybl+Npath-1,1)
    complex(kind(0d0)) Yn(Nsybl,1)
    complex(kind(0d0)) Y2(Nsybl,1)
    integer RdatI(Nsybl)
    integer RdatQ(Nsybl)
    complex(kind(0d0)) A(Nsybl,1)
    complex(kind(0d0)) Y2H(1,Nsybl)
    complex(kind(0d0)) R(1,1)
    complex(kind(0d0)) R2
    double precision EbN0(ConvSize)
    double precision BER(ConvSize)
    double precision Eig(1,Nsybl)
    double precision EbN0In(ConvSize)
    double precision AvRTNum(ConvSize)
    double precision AvUseChNum(ConvSize)
    double precision ConvStandard(ConvSize)

    !initialize
    Ampd=0.0d0
    Psig=0.0d0
    Pwgn=0.0d0
    Collect=0
    False=0
    loop=0
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
    RTNum=0
    UseChNum=0
    EbN0In=0.0d0
    AvRTNum=0.0d0
    AvUseChNum=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open (1, file='Cutoff.csv', status='replace')
    open (2, file='Cutoff_UseChNum.csv', status='replace')
    write(1,*) 'EbN0', ',', 'BER', ',', 'AvRTNum'

    !implimentation part
    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do
    
    !setup Convergence standard
    i=1
    do KEbN0=SEbN0, EEbN0, Step
        EbN0In(i) = 10.0**(dble(KEbN0)/10.0d0)
        if(KEbN0<0) then
            ConvStandard(i) = 1.0d0
        else
            ConvStandard(i) = 1.0d0*dexp(-0.23*dble(KEbN0))
        endif
        !print *, i, ConvStandard(i)
        i=i+1
    end do

    do loop=1, Nloop !Monte calro loop
        if(mod(loop,(Nloop)/10)==0) print *, loop

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
        call PPL(H,HE,Xppl,Eig,Nsybl,Npath,EbN0In,RTNum,UseChNum,ConvStandard,ConvSize,BERStandard,PV)

        do i=1, ConvSize
            !update AvRTNum
            AvRTNum(i) = AvRTNum(i) + dble(RTNum(i))/dble(Nloop)
            !update AvUseChNum
            AvUseChNum(i) = AvUseChNum(i) + dble(UseChNum(i))/dble(Nloop)
        end do

        !set information symbol
        S=0.0d0
        do i=1, Nsybl
            S(1,i) = cmplx(nint(rand())*2.0d0-1.0d0,nint(rand())*2.0d0-1.0d0,kind(0d0)) !-1or1
        end do
        TdatI=0
        TdatQ=0
        do i=1, Nsybl
            TdatI(i) = (real(S(1,i))+1)/2 !0or1
            TdatQ(i) = (aimag(S(1,i))+1)/2
        end do

        do k=1, ConvSize
            if(UseChNum(k)==0) cycle

            SU=(0.0d0,0.0d0)
            do i=1, UseChNum(k)
                do j=1, Nsybl
                    SU(j,i) = S(1,i) * PV(k,j,i)
                end do
            end do

            !transmit vector
            X = (0.0d0,0.0d0)
            do i=1,UseChNum(k)
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
            Noise = Noise / sqrt(2.0d0) * sqrt(1.0d0/EbN0In(k)/2.0d0)

            do i=1, Nsybl+Npath-1
                Psig(k) = Psig(k) + (abs(Y(i,1))**2.0d0)
                Pwgn(k) = Pwgn(k) + abs(Noise(i,1))**2.0d0
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
            do i=1, UseChNum(k)
                do j=1, Nsybl
                    A(j,1) = PV(k,j,i)
                end do
                !inner product
                call CAdjoint(Y2,Y2H,Nsybl,1)
                call CMultiply(Y2H,A,R,1,Nsybl,Nsybl,1)
                R(1,1) = conjg(R(1,1))
                R2 = R(1,1)
                if(real(R2)>0.0) then
                    RdatI(i) = 1
                elseif(real(R2)<0.0) then
                    RdatI(i) = 0
                endif
                if(aimag(R2)>0.0) then
                    RdatQ(i) = 1
                elseif(aimag(R2)<0.0) then
                    RdatQ(i) = 0
                endif
            end do

            !Bit Error Rate
            do i=1, UseChNum(k)
                if(RdatI(i)==TdatI(i)) then
                    Collect(k) = Collect(k) + 1
                elseif(RdatI(i)/=TdatI(i)) then
                    False(k) = False(k) + 1
                endif
                if(RdatQ(i)==TdatQ(i)) then
                    Collect(k) = Collect(k) + 1
                elseif(RdatQ(i)/=TdatQ(i)) then
                    False(k) = False(k) + 1
                endif
            end do
        end do
    end do
        
    do k=1, ConvSize
        if(AvUseChNum(k)==0.0d0) then
            EbN0(k) = 0.0d0
            BER(k) = 0.0d0
        else
            EbN0(k) = 10.0d0*dlog10(Psig(k)/Pwgn(k)/2.0d0)
            BER(k) = dble(False(k)) / (dble(Collect(k)) + dble(False(k)))
        endif
    end do

    do k=1, ConvSize
        write(2,*) nint(10.0d0*dlog10(EbN0In(k))), ',', AvUseChNum(k)
        write(1,*) EbN0(k), ',', BER(k), ',', AvRTNum(k)

        print *, 'EbN0=', EbN0(k)
        print *, 'BER=', BER(k)
        print *, 'AvRTNum=', AvRTNum(k)
        print *, 'UseChNum=', AvUseChNum(k)
        print *, ''
    end do

    !file close
    close(1)
    close(2)

    !time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program AvOth_BER