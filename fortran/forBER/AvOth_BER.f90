program AvOth_BER
    use PPLAvOth_BERmod
    use CALmod
    implicit none

    !-- run time declaration
    integer t1, t2, t_rate, t_max, diff

    !-- declaration
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=8
    integer,parameter :: KEbN0=10
    integer,parameter :: Nloop=10000
    double precision,parameter :: SAvOth=100.0d0
    double precision,parameter :: EAvOth=1.0d-5
    integer,parameter :: Step=10
    integer,parameter :: Vsize=8

    integer i,j,k
    double precision Threshold(1,Vsize)
    double precision Ampd(Npath,1)
    integer Collect(Vsize)
    integer False(Vsize)
    integer loop
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
    complex(kind(0d0)) Xppl(Nsybl,Nsybl)
    complex(kind(0d0)) PVO(Vsize,Nsybl,Nsybl)
    complex(kind(0d0)) HH(Nsybl,Nsybl+Npath-1)
    complex(kind(0d0)) HHH(Nsybl,Nsybl)
    complex(kind(0d0)) S(1,Nsybl)
    integer TdatI(1,Nsybl)
    integer TdatQ(1,Nsybl)
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
    double precision BER(Vsize)
    double precision PEO(Vsize,Nsybl)
    double precision ThrTMP

    !-- initialize
    Ampd(:,:)=0.0d0
    loop=0
    Cpath(:,:)=(0.0d0,0.0d0)
    H(:,:)=(0.0d0,0.0d0)
    HE(:,:)=(0.0d0,0.0d0)
    Xppl(:,:)=(0.0d0,0.0d0)
    HH(:,:)=(0.0d0,0.0d0)
    HHH(:,:)=(0.0d0,0.0d0)
    S(:,:)=0.0d0
    TdatI=0
    TdatQ=0
    SU(:,:)=(0.0d0,0.0d0)
    X(:,:)=(0.0d0,0.0d0)
    Pow=0.0d0
    Y(:,:)=(0.0d0,0.0d0)
    Noise(:,:)=(0.0d0,0.0d0)
    Yn(:,:)=(0.0d0,0.0d0)
    Y2(:,:)=(0.0d0,0.0d0)
    RdatI=0
    RdatQ=0
    A(:,:)=(0.0d0,0.0d0)
    Y2H(:,:)=(0.0d0,0.0d0)
    R(:,:)=(0.0d0,0.0d0)
    R2=(0.0d0,0.0d0)
    EbN0=0.0d0
    BER=0.0d0
    Threshold=10.0d0
    PVO=(0.0d0,0.0d0)
    PEO=0.0d0
    ThrTMP=0.0d0

    !-- time measurement start
    call system_clock(t1)

    !-- file open
    open (1, file='AvOth_BER(10dB)s32p8.csv', status='replace')

    !-- implimentation part
    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do

    !calculate EbN0
    EbN0 = 10.0d0**(dble(KEbN0)/10.0d0)

    ThrTMP = SAvOth
    do i=1, Vsize
        Threshold(1,i) = ThrTMP
        ThrTMP = ThrTMP *0.1d0
    end do

    Collect = 0
    False = 0
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
        call PPL(H,HE,Xppl,Nsybl,Npath,EbN0,Threshold,Vsize,PVO,PEO)

        do k=1, Vsize
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
                    SU(j,i) = S(1,i) * PVO(k,j,i)
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
                    A(j,1) = PVO(k,j,i)
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
                    Collect(k) = Collect(k) + 1
                elseif(RdatI(1,i)/=TdatI(1,i)) then
                    False(k) = False(k) + 1
                endif
                if(RdatQ(1,i)==TdatQ(1,i)) then
                    Collect(k) = Collect(k) + 1
                elseif(RdatQ(1,i)/=TdatQ(1,i)) then
                    False(k) = False(k) + 1
                endif
            end do
        end do
    end do
    
    do i=1, Vsize
        BER(i) = dble(False(i)) / dble(Collect(i)+False(i))
    end do

    do i=1, Vsize
        if(BER(i)>0.0) then
            write(1,*) Threshold(1,i), ',', BER(i)
            print *, 'Threshold=', Threshold(1,i)
            print *, 'BER=', BER(i)
        endif
    end do

    !file close
    close(1)

    !time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program AvOth_BER