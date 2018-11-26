program usePPLmod
    use PPLmod
    use CALmod
    implicit none

    !declaration
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=4
    integer,parameter :: PPLloop=1000

    integer i,j,k,l
    complex(kind(0d0)) X(Nsybl,Nsybl)
    complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
    complex(kind(0d0)) V(Nsybl,Nsybl)
    complex(kind(0d0)) NAISEKI(1,1)
    complex(kind(0d0)) V1(Nsybl)
    complex(kind(0d0)) V2(Nsybl)
    double precision Eig(1,Nsybl)
    double precision NAISEKI_TMP
    double precision AVGOTH

    !initialize
    X=(0.0d0,0.0d0)
    H=(0.0d0,0.0d0)
    HE=(0.0d0,0.0d0)
    NAISEKI=(0.0d0,0.0d0)
    NAISEKI_TMP=0.0d0
    AVGOTH=0.0d0
    Eig=0.0d0

    !-- file open
    open (1,file='usePPL.csv', status='replace')

    !set H
    do j=0, Npath-1
        do i=1, Nsybl
            H(i+j,i) = cmplx(0.1+0.1*j, 0.2+0.1*j, kind(0d0))
        end do
    end do

    !set HE
    do j=0, Npath-1
        do i=1, Nsybl+Npath-1
            HE(i+j,i) = cmplx(0.1+0.1*j, 0.2+0.1*j, kind(0d0))
        end do
    end do

	!set X
    do i=1, Nsybl
        do j=1, Nsybl
            X(i,j) = cmplx(1.0, 0.0, kind(0d0))
        end do
    end do

    do l=1, PPLloop
        call PPL(H,HE,X,Eig,Nsybl,Npath,1)

        !average othogonality of eiven vector
        NAISEKI(1,1) = cmplx(0.0,0.0,kind(0d0))
        do i=2, Nsybl
            !固有ベクトル群を１列のベクトルに格納
            do k=1, Nsybl
                V1(k) = X(k,i)
            end do
            do j=1, i-1
                !内積を取る固有ベクトルを格納
                do k=1, Nsybl
                    V2(k) = X(k,j)
                end do
                NAISEKI(1,1) = NAISEKI(1,1) + abs(dot_product(V1,V2))
            end do
        end do

        call CAbs(NAISEKI,NAISEKI_TMP,1,1)
        AVGOTH = NAISEKI_TMP

        if(mod(l,100)==0) then
            print *, l, AVGOTH
            write(1,*) l, ',', AVGOTH
        endif
    end do

    !-- file close
    close(1)
end program
