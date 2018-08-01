program eigenvaluePDF
    use CALmod
    use PPLmod
    implicit none

    !ifdef
    logical,parameter :: APPLY_PPL=.TRUE.

    !declatation
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=8
    integer,parameter :: PPLloop=500
    integer,parameter :: trial=10000

    double precision Ampd(Npath,1)
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
    complex(kind(0d0)) Xppl(Nsybl,Nsybl)
    complex(kind(0d0)) HH(Nsybl,Nsybl+Npath-1)
    complex(kind(0d0)) HHH(Nsybl,Nsybl)
    complex(kind(0d0)) V(Nsybl,Nsybl)
    double precision Eig(1,Nsybl)

    double precision,parameter :: stride=0.01d0
    double precision out
    double precision,allocatable :: result(:,:)
    integer i,j,loop
    double precision start,end,amp
    integer rank

    !initialization
    Ampd(:,:)=0.0d0
    Cpath(:,:)=(0.0d0,0.0d0)
    H(:,:)=(0.0d0,0.0d0)
    HE(:,:)=(0.0d0,0.0d0)
    Xppl(:,:)=(0.0d0,0.0d0)
    HH(:,:)=(0.0d0,0.0d0)
    HHH(:,:)=(0.0d0,0.0d0)
    V(:,:)=(0.0d0,0.0d0)
    Eig(:,:)=0.0d0

    out=0.0d0
    start=-5.0d0
    end=5.0d0
    amp=0.0d0
    rank=nint((end-start)/stride)+1

    !allocate initialize
    result(:,:)=0.0d0

    !file open
    open(1,file='fort_evPDF.csv', status='replace')

    !implementation
    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do



	do loop=1, trial
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

            call PPL(H,HE,Xppl,Eig,Nsybl,Npath,PPLloop)
        else
            call CSubstitute(V,HHH,Nsybl,Nsybl)
            call zdiag(Nsybl,V,Eig)
        endif

		out = normal()
		do j=1, rank
			if((out.ge.(start+stride*(j-1))).and.(out.lt.(start+stride*j))) then
				result(j,2) = result(j,2) + 1.0
			end if
		end do
	end do

    do i=1, rank
        result(i,1) = start + stride * (i-1)
        result(i,2) = result(i,2) / trial
    end do

    do i=1, rank
        write(1,*) result(i,1), ',', result(i,2)
    end do

    !file close
    close(1)

end program