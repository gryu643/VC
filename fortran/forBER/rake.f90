program rake
    use CALmod
    implicit none

    !-- run time declaration
    integer t1, t2, t_rate, t_max, diff

    !-- declaration
    integer,parameter :: SEbN0=1
    integer,parameter :: EEbN0=10
    integer,parameter :: Step=1
    integer,parameter :: Nloop=1000
    integer,parameter :: Nsybl=1
    integer,parameter :: Npath=2 !Npath>=2
    integer,parameter :: M_tapN=4

    integer i,j
    double precision Ampd(Npath,1)
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) H(2**M_tapN-1+Npath-1,2**M_tapN-1)
    double precision Ps
    double precision Pn
    integer Collect
    integer False
    integer loop
    integer KEbN0
    complex(kind(0d0)) S(Nsybl)
    complex(kind(0d0)) X(2**M_tapN-1,Nsybl)
    complex(kind(0d0)) Y(2**M_tapN-1+(11-1),Nsybl)
    complex(kind(0d0)) Y2(2**M_tapN-1+(11-1)+(2**M_tapN-1-1),Nsybl)
    complex(kind(0d0)) Y3(2**M_tapN-1+(11-1)+(2**M_tapN-1-1)+(11-1),Nsybl)
    integer TdatI(Nsybl)
    double precision Pow
    complex(kind(0d0)) Noise(2**M_tapN-1+(11-1),1)
    integer RdatI(Nsybl)
    double precision EbN0
    double precision BER
    integer M_weight(M_tapN)
    integer M_tap(M_tapN)
    complex(kind(0d0)) DSS_tap(2**M_tapN-1)
    complex(kind(0d0)) Rake_tap(Npath)
    double precision V(2**M_tapN-1)
    data M_weight/1,0,0,1/
    integer M_cal
    integer TMP
    complex(kind(0d0)) R
    complex(kind(0d0)) H_weight(11)
    complex(kind(0d0)) H_tap(11)

    !-- initialize
    Ps=0.0d0
    Pn=0.0d0
    Collect=0
    False=0
    loop=0
    KEbN0=0
    H=(0.0d0,0.0d0)
    S=(0.0d0,0.0d0)
    TdatI=0
    Pow=0.0d0
    Noise=(0.0d0,0.0d0)
    RdatI=0
    EbN0=0.0d0
    BER=0.0d0
    M_tap=0
    M_tap(1)=1
    M_cal=0
    R=(0.0d0,0.0d0)
    DSS_tap=(0.0d0,0.0d0)
    Rake_tap=(0.0d0,0.0d0)
    Y=(0.0d0,0.0d0)
    Y2=(0.0d0,0.0d0)
    Y3=(0.0d0,0.0d0)
    H_weight=(0.0d0,0.0d0)
    H_weight(1)=(1.0d0,0.0d0)
    H_weight(11)=(0.0d0,1.0d0)
    H_tap=(0.0d0,0.0d0)

    !-- time measurement start
    call system_clock(t1)

    !-- file open
    open (1, file='rake.csv', status='replace')
    open (2, file='rake_out.csv', status='replace')

    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
    end do

    !-- generate M-sequence
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
        V(i) =  TMP
        if(V(i)==0) V(i)=-1
        !print *, V(i)
    end do

    !-- implimentation part
    do KEbN0=SEbN0,EEbN0, Step !Ps/Pn loop
        Ps = 0.0d0
        Pn = 0.0d0
        Collect = 0
        False = 0

        do loop=1, Nloop !Monte calro loop
            !set H
!            do i=1, Npath
!                Cpath(i,1) = cmplx(normal(),normal(),kind(0d0))/sqrt(2.0d0)*Ampd(i,1)
!            end do
            Cpath(1,1)=cmplx(1.0d0,0.0d0,kind(0d0))
            Cpath(2,1)=cmplx(0.0d0,1.0d0,kind(0d0))

            !set H
            do i=1, 2**M_tapN-1
                do j=1, Npath
                    H(i+j-1,i) = Cpath(j,1)
                end do
            end do

            !set information symbol
            do i=1, Nsybl
                S(i) = cmplx(nint(rand())*2.0d0-1.0d0,0.0d0,kind(0d0)) !-1or1
            end do
            do i=1, Nsybl
                TdatI(i) = (real(S(i))+1)/2 !0or1
            end do

            !spreading
            do j=1, Nsybl
                do i=1, 2**M_tapN-1
                    X(i,j) = S(j) * V(i)
                end do
            end do

            !signal power = 1
            Pow = 0.0d0
            do i=1, 2**M_tapN-1
                Pow = Pow + abs(X(i,1))**2.0d0 / dble(2**M_tapN-1)
            end do
            X = X / sqrt(Pow)

            !generate Noise
            do i=1, 2**M_tapN-1+(11-1)
                Noise(i,1) = cmplx(normal(),normal(),kind(0d0))
            end do
            !apply EbN0
            Noise = Noise * sqrt(1.0d0/(10.0d0**(KEbN0/10.0d0))/2.0d0)

            !pass H --------------------------------------
            Y=(0.0d0,0.0d0)
            H_tap=(0.0d0,0.0d0)
            do i=1, 2**M_tapN-1+(11-1)
                !update despreading tap
                do j=11, 2, -1
                    H_tap(j) = H_tap(j-1)
                end do
                if(i>2**M_tapN-1) H_tap(1)=(0.0d0,0.0d0)
                if(i<=2**M_tapN-1) H_tap(1) = X(i,1)

                !calculate output
                do j=1, 11
                    Y(i,1) = Y(i,1) + H_tap(j)*H_weight(j)
                end do
            end do
            !---------------------------------------------

            !calculate Ps and Pn
            do i=1, 2**M_tapN-1+(11-1)
                Ps = Ps + abs(Y(i,1))**2.0d0
                Pn = Pn + abs(Noise(i,1))**2.0d0
            end do

            !add noise
            call CAdd(Y,Noise,Y,2**M_tapN-1+(11-1),1,2**M_tapN-1+(11-1),1)

            !despreading ----------------------------------
            Y2=(0.0d0,0.0d0)
            DSS_tap=(0.0d0,0.0d0)
            do i=1, 2**M_tapN-1+(11-1)+(2**M_tapN-1-1)
                !update despreading tap
                do j=2**M_tapN-1, 2, -1
                    DSS_tap(j) = DSS_tap(j-1)
                end do
                if(i>2**M_tapN-1+(11-1)) DSS_tap(1)=(0.0d0,0.0d0)
                if(i<=2**M_tapN-1+(11-1)) DSS_tap(1) = Y(i,1)

                !calculate output
                do j=1, 2**M_tapN-1
                    Y2(i,1) = Y2(i,1) + DSS_tap(j)*V(2**M_tapN-j)
                end do
            end do
            !----------------------------------------------

            !rake reception-----------------------------------------
            Y3=(0.0d0,0.0d0)
            H_tap=(0.0d0,0.0d0)
            do i=1, 2**M_tapN-1+(11-1)+(2**M_tapN-1-1)+(11-1)
                !update rake tap
                do j=11, 2, -1
                    H_tap(j) = H_tap(j-1)
                end do
                if(i>2**M_tapN-1+(11-1)+(2**M_tapN-1-1)) H_tap(1)=(0.0d0,0.0d0)
                if(i<=2**M_tapN-1+(11-1)+(2**M_tapN-1-1)) H_tap(1) = Y2(i,1)

                !calculate output
                do j=1, 11
                    Y3(i,1) = Y3(i,1) + H_tap(j)*conjg(H_weight(12-j))
                end do
            end do
            !----------------------------------------------

            !mabiki? --------------------------------------
            do i=1, 2**M_tapN-1+(11-1)+(2**M_tapN-1-1)
                !print *, i, Y3(i,1)
                write(2,*) i, ',', real(Y2(i,1))
                !print *, Y3(i,1)
            end do
            R = Y3(2**M_tapN-1+(11-1)+11,1)
            !----------------------------------------------

            RdatI = 0
            !Demodulation
            do i=1, Nsybl
                if(real(R)>0.0) then
                    RdatI(i) = 1
                elseif(real(R)<0.0) then
                    RdatI(i) = 0
                endif

                !Bit Error Rate
                if(RdatI(i)==TdatI(i)) then
                    Collect = Collect + 1
                elseif(RdatI(i)/=TdatI(i)) then
                    False = False + 1
                endif
            end do
        end do
        
        EbN0 = 10.0d0*dlog10(Ps/Pn/2.0d0)
        BER = dble(False) / (dble(Collect) + dble(False))
        if(BER>0.0) then
            write(1,*) KEbN0, ',', BER
            print *, 'EbN0=', EbN0
            print *, 'BER=', BER
        endif
    end do

    !-- file close
    close(1)
    close(2)

    !-- time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program