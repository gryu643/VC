program eigenvaluePDF_ideal
    implicit none

    !declaration
    double precision,parameter :: start=0.0d0
    double precision,parameter :: end=50.0d0
    double precision,parameter :: stride=0.01d0
    double precision lambda
    integer i
    integer rank
    double precision,allocatable :: result(:,:)



    !initializatioin
    rank=(nint((end-start)/stride))+1
    lambda=0.0d0

    !allocate
    allocate(result(rank,2))

    !file open
    open(1,file='evPDFideal.csv', status='replace')

    !implementation
    do i=1, rank
        lambda = start + stride*(i-1)
        result(i,1) = lambda
        result(i,2) = f(lambda)
    end do

    do i=1, rank
        write(1,*) result(i,1), ',', result(i,2)
    end do

    !file close
    close(1)

contains
    function f(lambda)
        implicit none
        double precision lambda
        double precision f

        f = (1.0d0/4.0d0)*((lambda**6)/36.0d0 &
        - (lambda**5)/2.0d0 + (lambda**4)*7.0d0/2.0d0 &
        - (lambda**3)*34.0d0/3.0d0 + (lambda**2)*18 - 12.0d0*lambda + 4.0d0)*exp(-lambda)

    end function

end program