program ch_est_err
    use ChEstmod
    implicit none

    !-- declaration
    integer,parameter :: M_tapN=4
    integer,parameter :: Npath=4
    integer,parameter :: EbN0=40
    integer,parameter :: Npilot=20
    integer,parameter :: Nloop=1000

    integer i,j,k
    integer loop
    integer pilot
    double precision Ampd(Npath,1)
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) H(2**M_tapN-1+Npath-1,2**M_tapN-1)
    complex(kind(0d0)) ChEst(Npath)
    double precision TMP
    double precision err

    !-- initialization
    pilot=0
    Ampd=0.0d0
    Cpath=(0.0d0,0.0d0)
    H=(0.0d0,0.0d0)
    ChEst=(0.0d0,0.0d0)
    TMP=0.0d0
    err=0.0d0

    !file open
    open (1, file='ch_est_err(EbN0=40).csv', status='replace')

    write(1,*) 'EbN0=', ',', EbN0

    !-- implementation
    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
    end do

    do pilot=1, Npilot
        err=0.0d0

        do loop=1, Nloop
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
            call ChEstimate(H,pilot,Npath,ChEst,M_tapN,EbN0)

            !calculate estimate err
            TMP=0.0d0
            do i=1, Npath
                TMP = TMP + abs(Cpath(i,1)-ChEst(i))/dble(Npath)
            end do
            err = err + TMP/dble(Nloop)
        end do

        !file output
        print *, "pilot=", pilot
        print *, "err=  ", err
        write(1,*) pilot, ',', err
    end do

    !file close
    close(1)

end program ch_est_err