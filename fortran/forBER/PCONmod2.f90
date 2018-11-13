module PCONmod2
    implicit none
contains
    subroutine Pcontrol2(lambda,EbN0,Pcon,Rep,Nsybl,info)
        use CALmod
        implicit none

        !--declaration
        integer Nsybl,i,Rep,info
        double precision lambda(1,Nsybl)
        double precision EbN0
        double precision Pcon(1,Nsybl)

        double precision mu
        double precision p1,p2
        double precision Min,Sum
        
        !--initialization
        mu=0.0d0
        p1=0.0d0
        p2=0.0d0
        info=1
        Min=0.0d0
        Sum=0.0d0

        !--implementation
        do i=1, Rep
            p1 = p1 + 1.0d0/lambda(1,i)
            p2 = p2 + 1.0d0/lambda(1,i)*dlog(dble(Rep)/lambda(1,i)/EbN0)
        end do

        mu = -EbN0/p1*(dble(Rep)+p2/EbN0)

        do i=1, Rep
            Pcon(1,i) = -1.0d0/lambda(1,i)/EbN0*(mu+dlog(dble(Rep)/lambda(1,i)/EbN0))
        end do

        do i=1, Rep
            if(0.0d0>Pcon(1,i).and.Pcon(1,i)<Min) then
                Min=Pcon(1,i)
                info=-1
            endif
        end do

        !set 0.01 for minimum value when minimum value is negative.
        if(info==-1) then
            do i=1, Rep
                Pcon(1,i) = Pcon(1,i) - (Min-0.01d0)
            end do
            !normalize Pcon
            do i=1, Rep
                Sum = Sum + Pcon(1,i)
            end do
            do i=1, Rep
                Pcon(1,i) = Pcon(1,i) / Sum * dble(Rep)
            end do
        endif
    end subroutine Pcontrol2
end module PCONmod2