module PCONmod3
    implicit none
contains
    subroutine Pcontrol3(lambda,EbN0,Pcon,Rep,Nsybl,info)
        use CALmod
        implicit none

        !--declaration
        integer Nsybl,i,Rep,info
        double precision lambda(1,Nsybl)
        double precision EbN0
        double precision Pcon(1,Nsybl)
        double precision A(Rep+1,Rep+1)
        double precision B(Rep+1,1)
        double precision C(Rep+1,1)
        double precision Min
        double precision Sum
        

        !--initialization
        A=0.0d0
        B=0.0d0
        C=0.0d0
        Min=0.0d0
        Sum=0.0d0

        !--implementation
        !AX=B
        !set A
        do i=1, Rep
            A(i,i) = lambda(1,i)*EbN0
            A(i,Rep+1) = 1.0d0
            A(Rep+1,i) = 1.0d0
        end do

        !set B
        do i=1, Rep
            B(i,1) = dlog(lambda(1,i)/dble(Rep)*EbN0)
        end do
        B(Rep+1,1) = dble(Nsybl)

        !calculate inverse matrix of A
        call InverseMat(A,Rep+1)

        !calculete Pcon
        call RMultiply(A,B,C,Rep+1,Rep+1,Rep+1,1)

        do i=1, Nsybl
            if(0.0d0>C(i,1)) then
                if(C(i,1)<Min) then
                    Min=C(i,1)
                    info=-1
                endif
            endif
        end do

        if(info==-1) then
            do i=1, Nsybl
                Pcon(i,1) = C(i,1) - (Min-0.01d0)
            end do
            do i=1, Nsybl
                Sum = Sum + Pcon(i,1)
            end do
            do i=1, Nsybl
                Pcon(i,1) = Pcon(i,1) / Sum * dble(Nsybl)
            end do
        else
            do i=1, Nsybl
                Pcon(i,1) = C(i,1)
            end do
        endif

    end subroutine Pcontrol3
end module PCONmod3