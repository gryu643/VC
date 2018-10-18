program VC_ideal_v2
    use PPLmod
    use CALmod
    implicit none

    !ifdef
    logical,parameter :: APPLY_PPL=.False.

    !run time declaration
    integer t1, t2, t_rate, t_max, diff

    !declaration
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=2
    integer,parameter :: SEbN0=-10
    integer,parameter :: EEbN0=40
    integer,parameter :: Step=10
    integer,parameter :: Nloop=100000
    integer,parameter :: PPLloop=500

    integer i,j
    integer loop
    integer KEbN0
    double precision Ampd(Npath,1)
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
    double precision LambdaEbN0(1,Nsybl)
    double precision AvLambdaEbN0
    double precision InstantBER(1,Nsybl)

    !initialize
    Ampd(:,:)=0.0d0
    loop=0
    KEbN0=0
    Cpath(:,:)=(0.0d0,0.0d0)
    H(:,:)=(0.0d0,0.0d0)
    HE(:,:)=(0.0d0,0.0d0)
    Xppl(:,:)=(0.0d0,0.0d0)
    HH(:,:)=(0.0d0,0.0d0)
    HHH(:,:)=(0.0d0,0.0d0)
    V(:,:)=(0.0d0,0.0d0)
    EbN0=0.0d0
    BER=0.0d0
    Eig=0.0d0
    LambdaEbN0=0.0d0
    InstantBER=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open (1, file='VC_ideal_v2(s32p2).csv', status='replace')

    !implimentation part
    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i,1) = sqrt(1.0d0/(2.0d0**(dble(i)-1.0d0))) !Exp. atten.
    end do

    do KEbN0=SEbN0, EEbN0, Step !Eb/N0 loop
        EbN0 = 10.0d0**(dble(KEbN0)/10.0d0)
        BER=0.0d0
        AvLambdaEbN0=0.0d0

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
                if(Npath==1) then
                    do i=1, Nsybl
                        Eig(1,i) = HHH(i,i)
                        V(i,i) = cmplx(1.0d0,0.0d0,kind(0d0))
                    end do
                else
                    call CSubstitute(V,HHH,Nsybl,Nsybl)
                    call decomp_zheevd(Nsybl,V,Eig)
    !                call decomp_zheev(Nsybl,V,Eig)
    !                call decomp_zgeev(Nsybl,V,Eig)
    !                call decomp_zhpev(Nsybl,V,Eig)
                endif
            endif

            call sort(Eig,Nsybl)

            do i=1, Nsybl
                LambdaEbN0(1,i) = Eig(1,i)*EbN0
                InstantBER(1,i) = 1.0d0/2.0d0*erfc(sqrt(LambdaEbN0(1,i)))
                BER = BER + InstantBER(1,i)/dble(Nsybl)/dble(Nloop)
                AvLambdaEbN0 = AvLambdaEbN0 + LambdaEbN0(1,i)/dble(Nsybl)/dble(Nloop)
            end do
        end do

        EbN0 = 10.0d0*dlog10(AvLambdaEbN0) !QPSK rate =2
        select case(Npath)
            case(1)
                EbN0 = EbN0+0.091891291
            case(2)
                EbN0 = EbN0-0.213893347
            case(4)
                EbN0 = EbN0-0.328150749
            case(8)
                EbN0 = EbN0-0.893632232
        end select

        if(BER>0.0) then
            write(1,*) EbN0, ',', BER
            print *, 'EbN0=', EbN0
            print *, 'BER=', BER
        endif
    end do

    !file close
    close(1)

    !time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program VC_ideal_v2