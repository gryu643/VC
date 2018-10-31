module PCONmod
    implicit none
contains
    subroutine Pcontroll(lambda,EbN0,Pcon,Nsymbl)
        use CALmod
        implicit none

        !--declaration
        integer Nsymbl,i
        double precision lambda(1,Nsymbl)
        double precision EbN0
        double precision Pcon(1,Nsymbl)
        double precision A(Nsymbl+1,Nsymbl+1)
        double precision B(Nsymbl+1,1)
        double precision C(Nsymbl+1,1)

        !--initialization
        A=0.0d0
        B=0.0d0
        C=0.0d0

        !--implementation
        !AX=B
        !set A
        do i=1, Nsymbl
            A(i,i) = lambda(1,i)*EbN0
            A(i,Nsymbl+1) = 1.0d0
            A(Nsymbl+1,i) = 1.0d0
        end do

        !set B
        
        do i=1, Nsymbl
            B(i,1) = dlog(lambda(1,i)/Nsymbl*EbN0)
        end do
        B(Nsymbl+1,1) = 1.0d0

        !calculate inverse matrix of A
        call InverseMat(A,Nsymbl+1)

        !calculete Pcon
        call RMultiply(A,B,C,Nsymbl+1,Nsymbl+1,Nsymbl+1,1)

        do i=1, Nsymbl
            Pcon(1,i) = C(i,1)
        end do
    end subroutine Pcontroll
end module PCONmod