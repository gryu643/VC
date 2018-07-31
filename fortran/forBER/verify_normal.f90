program verify_normal
    use CALmod
    implicit none

    !declaration
    integer,parameter :: trial=1000000
    double precision,parameter :: stride=0.01d0
    double precision out
    double precision,allocatable :: result(:,:)
    integer i,j
    double precision sum,avg,var
    double precision start,end,amp
    integer rank

    !initialization
    out=0.0d0
    sum=0.0d0
    avg=0.0d0
    var=0.0d0
    start=-5.0d0
    end=5.0d0
    amp=0.0d0
    rank=nint((end-start)/stride)+1
    
    !allocate
    allocate(result(rank,2))

    !allocate initialize
    result(:,:)=0.0d0

    !file open
    open(1,file='fortran_normal.csv', status='replace')


    !implimentation
	do i=1, trial
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

    !calculate avg
    do i=1,rank
        avg = avg + result(i,1)*result(i,2)
    end do

    !calculate var
    do i=1, rank
        var = var + ((result(i,1)-avg)**2)*result(i,2)
    end do

    print *, 'avg...', avg
    print *, 'var...', var
    !file close
    close(1)
end program