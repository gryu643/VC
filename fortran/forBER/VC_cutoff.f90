program VC_cutoff
    use PPLmod
    use CALmod
    implicit none

    !ifdef
    logical,parameter :: APPLY_PPL=.False.

    !run time declaration
    integer t1, t2, t_rate, t_max, diff

    !declaration
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=8
    integer,parameter :: SEbN0=-10
    integer,parameter :: EEbN0=40
    integer,parameter :: Step=10
    integer,parameter :: Nloop=10000
    integer,parameter :: PPLloop=500

    integer i,j
    double precision Ampd(Npath,1)
    double precision Psig(1,Nsybl)
    double precision Pwgn(1,Nsybl)
    integer Collect(1,Nsybl)
    integer False(1,Nsybl)
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
    double precision EbN0(1,Nsybl)
    double precision BER(1,Nsybl)
    double precision Eig(1,Nsybl)
    integer KRep
    complex(kind(0d0)) TMP(Nsybl)

    !initialize
    Ampd(:,:)=0.0d0
    Psig=0.0d0
    Pwgn=0.0d0
    Collect=0
    False=0
    loop=0
    KEbN0=0
    Cpath(:,:)=(0.0d0,0.0d0)
    H(:,:)=(0.0d0,0.0d0)
    HE(:,:)=(0.0d0,0.0d0)
    Xppl(:,:)=(0.0d0,0.0d0)
    HH(:,:)=(0.0d0,0.0d0)
    HHH(:,:)=(0.0d0,0.0d0)
    S(:,:)=0.0d0
    TdatI=0
    TdatQ=0
    V(:,:)=(0.0d0,0.0d0)
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
    Eig=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open (7, file='VCcutoff(s32p8)ch1.csv', status='replace')
    open (8, file='VCcutoff(s32p8)ch2.csv', status='replace')
    open (9, file='VCcutoff(s32p8)ch3.csv', status='replace')
    open (10, file='VCcutoff(s32p8)ch4.csv', status='replace')
    open (11, file='VCcutoff(s32p8)ch5.csv', status='replace')
    open (12, file='VCcutoff(s32p8)ch6.csv', status='replace')
    open (13, file='VCcutoff(s32p8)ch7.csv', status='replace')
    open (14, file='VCcutoff(s32p8)ch8.csv', status='replace')
    open (15, file='VCcutoff(s32p8)ch9.csv', status='replace')
    open (16, file='VCcutoff(s32p8)ch10.csv', status='replace')
    open (17, file='VCcutoff(s32p8)ch11.csv', status='replace')
    open (18, file='VCcutoff(s32p8)ch12.csv', status='replace')
    open (19, file='VCcutoff(s32p8)ch13.csv', status='replace')
    open (20, file='VCcutoff(s32p8)ch14.csv', status='replace')
    open (21, file='VCcutoff(s32p8)ch15.csv', status='replace')
    open (22, file='VCcutoff(s32p8)ch16.csv', status='replace')
    open (23, file='VCcutoff(s32p8)ch17.csv', status='replace')
    open (24, file='VCcutoff(s32p8)ch18.csv', status='replace')
    open (25, file='VCcutoff(s32p8)ch19.csv', status='replace')
    open (26, file='VCcutoff(s32p8)ch20.csv', status='replace')
    open (27, file='VCcutoff(s32p8)ch21.csv', status='replace')
    open (28, file='VCcutoff(s32p8)ch22.csv', status='replace')
    open (29, file='VCcutoff(s32p8)ch23.csv', status='replace')
    open (30, file='VCcutoff(s32p8)ch24.csv', status='replace')
    open (31, file='VCcutoff(s32p8)ch25.csv', status='replace')
    open (32, file='VCcutoff(s32p8)ch26.csv', status='replace')
    open (33, file='VCcutoff(s32p8)ch27.csv', status='replace')
    open (34, file='VCcutoff(s32p8)ch28.csv', status='replace')
    open (35, file='VCcutoff(s32p8)ch29.csv', status='replace')
    open (36, file='VCcutoff(s32p8)ch30.csv', status='replace')
    open (37, file='VCcutoff(s32p8)ch31.csv', status='replace')
    open (38, file='VCcutoff(s32p8)ch32.csv', status='replace')
    !implimentation pa
    !channel gain paramet
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do

    do KEbN0=SEbN0, EEbN0, Step !Eb/N0 loop
        Psig = 0.0d0
        Pwgn = 0.0d0
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

            if(APPLY_PPL) then
                !set Xppl
                do i=1, Nsybl
                    do j=1, Nsybl
                        Xppl(i,j) = cmplx(1.0d0, 0.0d0, kind(0d0))
                    end do
                end do
            endif

            !set HH
            call CAdjoint(H,HH,Nsybl+Npath-1,Nsybl)

            !set HHH
            call CMultiply(HH,H,HHH,Nsybl,Nsybl+Npath-1,Nsybl+Npath-1,Nsybl)

            !eigenvalue decomposition
            if(APPLY_PPL) then
                call PPL(H,HE,Xppl,Eig,Nsybl,Npath,PPLloop)
                call CSubstitute(V,Xppl,Nsybl,Nsybl)
            else
                call CSubstitute(V,HHH,Nsybl,Nsybl)
                call decomp_zheevd(Nsybl,V,Eig)
!                call decomp_zheev(Nsybl,V,Eig)
!                call decomp_zgeev(Nsybl,V,Eig)
!                call decomp_zhpev(Nsybl,V,Eig)

                !ライブラリ出力時は逆順なので、ソート
                call sort(Eig,Nsybl)
                do i=1, Nsybl/2
                    do j=1, Nsybl
                        TMP(j) = V(j,i)
                        V(j,i) = V(j,Nsybl+1-i)
                        V(j,Nsybl+1-i) = TMP(j)
                    end do
                end do
            endif

            do KRep=1, Nsybl

                !set information symbol
                do i=1, KRep
                    S(1,i) = cmplx(nint(rand())*2.0d0-1.0d0,0.0d0,kind(0d0)) !-1or1
                end do
                do i=1, KRep
                    TdatI(1,i) = (real(S(1,i))+1)/2 !0or1
                    TdatQ(1,i) = (aimag(S(1,i))+1)/2
                end do
                do i=1, KRep
                    do j=1, Nsybl
                        SU(j,i) = S(1,i) * V(j,i)
                    end do
                end do

                !transmit vector
                X = (0.0d0,0.0d0)
                do i=1,KRep
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
                    Psig(1,Krep) = Psig(1,KRep) + (abs(Y(i,1))**2.0d0)
                    Pwgn(1,KRep) = Pwgn(1,KRep) + abs(Noise(i,1))**2.0d0
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
                        Collect(1,KRep) = Collect(1,KRep) + 1
                    elseif(RdatI(1,i)/=TdatI(1,i)) then
                        False(1,KRep) = False(1,KRep) + 1
                    endif
                    if(RdatQ(1,i)==TdatQ(1,i)) then
                        Collect(1,KRep) = Collect(1,KRep) + 1
                    elseif(RdatQ(1,i)/=TdatQ(1,i)) then
                        False(1,KRep) = False(1,KRep) + 1
                    endif
                end do
            end do
        end do
        
        do i=1,Nsybl
            EbN0(1,i) = 10.0d0*dlog10(Psig(1,i)/Pwgn(1,i)) !QPSK rate = 1
            BER(1,i) = dble(False(1,Nsybl)) / (dble(Collect(1,Nsybl)) + dble(False(1,Nsybl)))
        end do
        do i=1,Nsybl
            if(BER(1,i)>0.0) then
                write(i+6,*) EbN0(1,i), ',', BER(1,i)
                print *, 'Symbl=', i, 'EbN0=', EbN0(1,i)
                print *, 'Symbl=', i, 'BER=', BER(1,i)
            endif
        end do
    end do

    !file close
    do i=7, 39
        close(i)
    end do

    !time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program