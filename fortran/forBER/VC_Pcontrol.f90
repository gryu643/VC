program VC_Pcontrol
    use PPLmod
    use CALmod
    use PCONmod
    use PCONmod2
    implicit none

    !ifdef
    logical,parameter :: APPLY_PPL=.False.

    !run time declaration
    integer t1, t2, t_rate, t_max, diff

    !declaration
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=4
    integer,parameter :: SEbN0=-10
    integer,parameter :: EEbN0=40
    integer,parameter :: Step=10
    integer,parameter :: Nloop=10000
    integer,parameter :: PPLloop=500

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
    double precision Pt(1,Nsybl)
    complex(kind(0d0)) TMP
    double precision EbN0t
    double precision BERMin
    double precision UsePt(1,Nsybl)
    integer info
    integer TNum
    double precision AvTNum
    double precision Sum
    double precision BER2
    double precision AvBER
    double precision AvBER2
    double precision AvEbN0
    double precision AvEbN02

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
    Pt=0.0d0
    TMP=0.0d0
    UsePt=0.0d0
    info=1

    !time measurement start
    call system_clock(t1)

    !file open
    open (1, file='VC_Pcon(s32p4).csv', status='replace')
    open (2, file='VC_Pcon(s32p4)ideal.csv', status='replace')

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
        AvTNum=0.0d0
        AvBER=0.0d0
        AvBER2=0.0d0
        AvEbN0=0.0d0
        AvEbN02=0.0d0

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
            endif

            !sort Eig descending order
            call sort(Eig,Nsybl)
            do i=1,Nsybl
                do j=1, Nsybl
                    TMP = V(j,i)
                    V(j,i) = V(j,Nsybl+1-i)
                    V(j,Nsybl+1-i) = TMP
                end do
            end do
            !transmit power control
            EbN0t = 10**(KEbN0/10.0d0)
            do i=Nsybl, 1, -1
                Pt=0.0d0
                info=1
                call Pcontrol(Eig,EbN0t,Pt,i,Nsybl,info)
                if(info==1) then
                    do j=1, Nsybl
                        UsePt(1,j) = Pt(1,j)
                    end do
                    TNum=i
                    exit
                endif
            end do
            AvTNum = AvTNum + dble(TNum)/dble(Nloop)
            if(info==-1) cycle
            !-- Pcontrol2
!            UsePt=0.0d0
!            info = 1
!            call Pcontrol2(Eig,EbN0t,UsePt,Nsybl,Nsybl,info)
!            if(info==-1) cycle
!            TNum = info

            !set information symbol
            S=(0.0d0,0.0d0)
            do i=1, TNum
                S(1,i) = cmplx(nint(rand())*2.0d0-1.0d0,nint(rand())*2.0d0-1.0d0,kind(0d0)) !-1or1
            end do

            TdatI=0
            TdatQ=0
            do i=1, TNum
                TdatI(1,i) = (real(S(1,i))+1)/2 !0or1
                TdatQ(1,i) = (aimag(S(1,i))+1)/2
            end do
            SU=(0.0d0,0.0d0)
            do i=1, TNum
                do j=1, Nsybl
                    SU(j,i) = S(1,i) * V(j,i) * sqrt(UsePt(1,i))
                end do
            end do

            do i=1, Nsybl
                AvEbN0 = AvEbN0 + EbN0t*Eig(1,i)/dble(Nsybl)/dble(Nloop)
            end do

            do i=1, TNum
                AvEbN02 = AvEbN02 + EbN0t*Eig(1,i)*UsePt(1,i)/dble(TNum)/dble(Nloop)
            end do

            BER=0.0d0
            do i=1, Nsybl
                BER=BER+1.0d0/2.0d0*erfc(sqrt(EbN0t*Eig(1,i)))/dble(Nsybl)
            end do
            AvBER = AvBER + BER/dble(Nloop)

            BER2=0.0d0
            do i=1, TNum
                BER2=BER2+1.0d0/2.0d0*erfc(sqrt(EbN0t*Eig(1,i)*UsePt(1,i)))/dble(TNum)
            end do
            AvBER2 = AvBER2 + BER2/dble(Nloop)

            if(BER>=BER2) then
                !print *, 'KEbN0, BER1, BER2=', KEbN0, BER, BER2
            else
                !print *, 'KEbN0, BER1, BER2=', KEbN0, BER, BER2, 'BER1<BER2'
                do i=1, TNum
                    !print *, i, UsePt(1,i)
                end do
            endif

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
                Pow = Pow + (abs(X(i,1))**2.0d0)/dble(Nsybl)
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
            do i=1, TNum
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
        AvEbN0 = 10.0d0*dlog10(AvEbN0)
        AvEbN02 = 10.0d0*dlog10(AvEbN02)
        if(BER>0.0) then
            write(1,*) AvEbN0, ',', AvBER, ',', AvEbN02, ',', AvBER2, ',', AvTNum
            print *, 'EbN0=', EbN0
            print *, 'BER=', BER
            print *, 'IdealBER=', AvBER
        endif
    end do

    !file close
    close(1)
    close(2)

    !time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program