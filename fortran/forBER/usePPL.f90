program usePPL
    use PPLmod
    use PPLNmod
    use CALmod
    use ChEstmod
    implicit none

    !declaration
    logical,parameter :: ChEstimate_=.True.
    logical,parameter :: PPLNoise_=.True.

    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=4
    integer,parameter :: PPLloop=1000
    integer,parameter :: M_tapN=4
    integer,parameter :: Npilot=4
    integer,parameter :: EbN0=10

    integer i,j,k,l
    complex(kind(0d0)) H_forChEst(2**M_tapN-1+Npath-1,2**M_tapN-1)
    double precision Ampd(Npath,1)
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) X(Nsybl,Nsybl)
    complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
    complex(kind(0d0)) V(Nsybl,Nsybl)
    double precision Eig(1,Nsybl)
    complex(kind(0d0)) ChEst(Npath)
    complex(kind(0d0)) H_est(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HH(Nsybl,Nsybl+Npath-1)
    complex(kind(0d0)) HHH(Nsybl,Nsybl)
    complex(kind(0d0)) TMP
    double precision AvOth
    double precision Pow(Nsybl)

    !initialize
    Ampd=0.0d0
    X=(0.0d0,0.0d0)
    H=(0.0d0,0.0d0)
    HH=(0.0d0,0.0d0)
    HHH=(0.0d0,0.0d0)
    HE=(0.0d0,0.0d0)
    Eig=0.0d0
    Cpath=(0.0d0,0.0d0)
    ChEst=(0.0d0,0.0d0)
    H_est=(0.0d0,0.0d0)
    TMP=(0.0d0,0.0d0)
    AvOth=0.0d0

    !-- file open
    open (1,file='usePPL.csv', status='replace')

    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
    end do

    !set Cpath
    do i=1, Npath
        Cpath(i,1) = cmplx(normal(),normal(),kind(0d0))/sqrt(2.0d0)*Ampd(i,1)
    end do

    !set H
    H=(0.0d0,0.0d0)
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

    !set H for channel estimate
    do i=1, 2**M_tapN-1
        do j=1, Npath
            H_forChEst(i+j-1,i) = Cpath(j,1)
        end do
    end do

    if(ChEstimate_) then
        !channel estimate
        call ChEstimate(H_forChEst,Npilot,Npath,ChEst,M_tapN,EbN0)

        !make H_est from CSI
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
        call CSubstitute(X,HHH,Nsybl,Nsybl)
        call decomp_zheevd(Nsybl,X,Eig)
        !------------------------------------------------------------

        !sort
        call sort(Eig,Nsybl)
        do i=1,Nsybl/2
            do j=1, Nsybl
                TMP = X(j,i)
                X(j,i) = X(j,Nsybl+1-i)
                X(j,Nsybl+1-i) = TMP
            end do
        end do

        Pow = 0.0d0
        do j=1, Nsybl
            do i=1, Nsybl
                Pow(j) = Pow(j) + (abs(X(i,j))**2.0d0)/Nsybl
            end do
            do i=1, Nsybl
                X(i,j) = X(i,j) / sqrt(Pow(j))
            end do
        end do
    else
        !set X
        do i=1, Nsybl
            do j=1, Nsybl
                X(i,j) = cmplx(1.0, 0.0, kind(0d0))
            end do
        end do
    endif

    do l=1, PPLloop
        !call PPL(H,HE,X,Eig,Nsybl,Npath,1)
        call PPLN(H,HE,X,Eig,Nsybl,Npath,EbN0,1,PPLNoise_)

        if(mod(l,(PPLloop/10))==0) then
            AvOth = orthogonal(X,Nsybl)
            print *, l, AvOth
            write(1,*) l, ',', AvOth
        endif
    end do

    !-- file close
    close(1)
end program usePPL
