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
        double precision parts1,parts2
        integer flag
        

        !--initialization
        mu=0.0d0
        parts1=0.0d0
        parts2=0.0d0
        flag=0

        !--implementation
        do i=1, Nsybl
            parts1 = parts2 + 1.0d0/lambda(1,i)/EbN0*dlog(lambda(1,i)/dble(Nsybl)*EbN0)
        end do
        do i=1, Nsybl
            parts2 = parts2 + 1.0d0/lambda(1,i)/EbN0
        end do

        mu = -1.0d0*dexp((parts1-1.0d0)/parts2)
        do i=1, Nsybl
            if(lambda(1,i)/dble(Nsybl)*EbN0>=-mu) then
                Pcon(1,i) = 1.0d0/lambda(1,i)/EbN0*(dlog(lambda(1,i)/dble(Nsybl)*EbN0)-dlog(-mu))
            else
                Pcon(1,i) = 0.0d0
                flag=flag+1
            endif
        end do
        if(flag==Nsybl) then
            info=-1
        else
            info=Nsybl-flag
        endif
    end subroutine Pcontrol2
end module PCONmod2