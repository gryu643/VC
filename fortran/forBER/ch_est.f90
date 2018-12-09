program ch_est
    use CALmod
    use ChEstmod
    implicit none

    !-- run time declaration
    integer t1, t2, t_rate, t_max, diff

    !-- declaration
    integer,parameter :: SEbN0=-10
    integer,parameter :: EEbN0=40
    integer,parameter :: Step=10
    integer,parameter :: Nloop=10000
    integer,parameter :: Nsybl=32
    integer,parameter :: Npilot=16
    integer,parameter :: Npath=4 !Npath>=2
    integer,parameter :: M_tapN=4

    integer i,j
    complex(kind(0d0)) ChEst(Npath)
    double precision Psig
    double precision Pwgn
    double precision Ampd(Npath,1)
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) H(2**M_tapN-1+Npath-1,2**M_tapN-1)
    integer Collect
    integer False
    integer loop
    integer KEbN0
    double precision Pow
    double precision EbN0
    double precision BER
    complex(kind(0d0)) H_est(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HH(Nsybl,Nsybl+Npath-1)
    complex(kind(0d0)) HHH(Nsybl,Nsybl)
    complex(kind(0d0)) S(1,Nsybl)
    integer TdatI(1,Nsybl)
    integer TdatQ(1,Nsybl)
    complex(kind(0d0)) V(Nsybl,Nsybl)
    complex(kind(0d0)) SU(Nsybl,Nsybl)
    complex(kind(0d0)) X(Nsybl,Nsybl)
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
    double precision Eig(1,Nsybl)

    !-- initialize
    Ampd=0.0d0
    Psig=0.0d0
    Pwgn=0.0d0
    Collect=0
    False=0
    loop=0
    KEbN0=0
    Cpath=(0.0d0,0.0d0)
    H=(0.0d0,0.0d0)
    H_est=(0.0d0,0.0d0)
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

    !-- time measurement start
    call system_clock(t1)

    !-- file open
    open (1, file='ch_est(Npilot=16).csv', status='replace')

    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
    end do

    !-- implimentation part
    do KEbN0=SEbN0,EEbN0, Step !Ps/Pn loop
        Psig = 0.0d0
        Pwgn = 0.0d0
        Collect = 0
        False = 0

        do loop=1, Nloop !Monte calro loop
            !set H
            do i=1, Npath
                Cpath(i,1) = cmplx(normal(),normal(),kind(0d0))/sqrt(2.0d0)*Ampd(i,1)
            end do

            !set H
            do i=1, 2**M_tapN-1
                do j=1, Npath
                    H(i+j-1,i) = Cpath(j,1)
                end do
            end do

            !channel estimation
            call ChEstimate(H,Npilot,Npath,ChEst,M_tapN,KEbN0)

            !make H from CSI
            do i=1, Nsybl
                do j=1, Npath
                    H_est(i+j-1,i) = Chest(j) !Cpath(j,1)
                end do
            end do

            !set HH
            call CAdjoint(H_est,HH,Nsybl+Npath-1,Nsybl)

            !set HHH
            call CMultiply(HH,H_est,HHH,Nsybl,Nsybl+Npath-1,Nsybl+Npath-1,Nsybl)

            !QP decomposition -------------------------------------------
            call CSubstitute(V,HHH,Nsybl,Nsybl)
            call decomp_zheevd(Nsybl,V,Eig)
            !------------------------------------------------------------

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

            call CMultiply(H_est,X,Y,Nsybl+Npath-1,Nsybl,Nsybl,1) !pass H

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
       end do
        
        EbN0 = 10.0d0*dlog10(Psig/Pwgn/2.0d0)
        BER = dble(False) / (dble(Collect) + dble(False))
       if(BER>0.0) then
            write(1,*) EbN0, ',', BER
            print *, 'EbN0=', EbN0
            print *, 'BER=', BER
        endif
    end do

    !-- file close
    close(1)

    !-- time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program ch_est