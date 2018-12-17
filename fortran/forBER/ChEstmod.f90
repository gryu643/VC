module ChEstmod
    !-- Channel Estimation --!
    use CALmod
    implicit none
contains
    subroutine ChEstimate(H,Np,Npath,ChEst,M_tapN,EbN0)
        implicit none

        !-- argument
        integer Np
        integer Npath
        integer M_tapN
        integer EbN0
        complex(kind(0d0)) ChEst(Npath)
        complex(kind(0d0)) H(2**M_tapN-1+Npath-1,2**M_tapN-1)
        
        !-- declaration
        integer i,j,k
        integer M_weight(M_tapN)
        integer M_tap(M_tapN)
        complex(kind(0d0)) DSS_tap(2**M_tapN-1)
        double precision V(2**M_tapN-1)
        integer M_cal
        integer TMP
        complex(kind(0d0)) S(Np)
        complex(kind(0d0)) X(2**M_tapN-1,Np)
        complex(kind(0d0)) Y(2**M_tapN-1+Npath-1,Np)
        complex(kind(0d0)) Y2(2**M_tapN-1+Npath-1,Np)
        complex(kind(0d0)) Noise(2**M_tapN-1+Npath-1,Np)

        !-- initialization
        M_tap=0
        M_tap(1)=1 !set initial value, [1 0 0 0]
        M_cal=0
        DSS_tap=(0.0d0,0.0d0)
        V=0.0d0
        TMP=0
        S=(0.0d0,0.0d0)
        X=(0.0d0,0.0d0)
        Y=(0.0d0,0.0d0)
        Y2=(0.0d0,0.0d0)
        Noise=(0.0d0,0.0d0)
        M_weight=(/1,0,0,1/)

        !-- implementation
        !generate M-sequence
        do i=1, 2**M_tapN-1
            !calculate xor
            M_cal = M_tap(M_tapN)*M_weight(M_tapN)
            do j=M_tapN-1, 1, -1
                M_cal = abs(M_cal-M_tap(j)*M_weight(j))
            end do

            !update M_tap
            do j=1, M_tapN
                TMP = M_tap(j)
                M_tap(j) = M_cal
                M_cal = TMP
            end do

            !output M-sequence
            V(i) =  TMP
            if(V(i)==0) V(i)=-1
            !print *, V(i)
        end do

        !set pilot symbol
        do i=1, Np
            S(i) = cmplx(1.0d0,0.0d0,kind(0d0)) !-1or1
        end do

        !spreading
        do j=1, Np
            do i=1, 2**M_tapN-1
                X(i,j) = S(j) * V(i)
            end do
        end do

        !generate Noise
        do j=1, Np
            do i=1, 2**M_tapN-1+Npath-1
                Noise(i,j) = cmplx(normal(),normal(),kind(0d0))
            end do
        end do
        !apply EbN0
        Noise = Noise / sqrt(2.0d0) * sqrt(1.0d0/(10.0d0**(EbN0/10.0d0))/2.0d0)

        !pass H --------------------------------------
        call CMultiply(H,X,Y,2**M_tapN-1+Npath-1,2**M_tapN-1,2**M_tapN-1,Np)
        !---------------------------------------------

        !add noise
        call CAdd(Y,Noise,Y,2**M_tapN-1+Npath-1,1,2**M_tapN-1+Npath-1,1)
        
        !add pilot(to decrease noise)
        if(Np>1) then
            do i=2, Np
                do j=1, 2**M_tapN-1+Npath-1
                    Y(j,1) = Y(j,1) + Y(j,i)
                end do
            end do
            do j=1, 2**M_tapN-1+Npath-1
                Y(j,1) = Y(j,1) / dble(Np)
            end do
        endif

        !despreading ----------------------------------
        Y2=(0.0d0,0.0d0)
        DSS_tap=(0.0d0,0.0d0)
        do i=1, 2**M_tapN-1+Npath-1
            !update despreading tap
            do j=2**M_tapN-1, 2, -1
                DSS_tap(j) = DSS_tap(j-1)
            end do
            DSS_tap(1) = Y(i,1)

            !calculate output
            do j=1, 2**M_tapN-1
                Y2(i,1) = Y2(i,1) + DSS_tap(j)*V(2**M_tapN-j)
            end do
        end do
        !----------------------------------------------

        !channel estimation ---------------------------
        do i=1, Npath
            ChEst(i) = Y2(2**M_tapN-1-1+i,1) / dble(2**M_tapN-1)
        end do
        !----------------------------------------------

    end subroutine ChEstimate
end module ChEstmod