program usePPLmod
    use usePPLmod
    implicit none

    integer :: Nsybl=16
    integer :: Npath=8
    integer :: PPLloop=500

    complex(kind(0d0)) :: X(Nsybl,Nsybl)=(0.0,0.0)
    complex(kind(0d0)) :: H(Nsybl+Npath-1,Nsybl)=(0.0,0.0)
    complex(kind(0d0)) :: HE(Nsybl+2*(Npath-1)),Nsybl+Npath-1)=(0.0,0.0)
    complex(kind(0d0)) :: V(Nsybl,Nsybl)=(0.0,0.0)

    complex(kind(0d0)) :: NAISEKI(1,1)=(0.0,0.0)
    complex(kind(0d0)),allocatable :: U(:,:)
    complex(kind(0d0)),allocatable :: N(:,:)
    complex(kind(0d0)),allocatable :: UH(:,:)
    complex(kind(0d0)) :: UHN(1,1)=(0.0)
    double precision :: NAISEKI_TMP=0.0
    double precision :: sC2=0.0
    double precision :: AVGOTH=0.0

    !set H
    do j=0, PATH-1
        do i=1, SYMBL
            H(i+j,i) = cmplx(0.1+0.1*j, 0.2+0.1*j, kind(0d0))
        end do
    end do

    !set HE
    do j=0, PATH-1
        do i=1, SYMBL+PATH-1
            HE(i+j,i) = cmplx(0.1+0.1*j, 0.2+0.1*j, kind(0d0))
        end do
    end do

	!set X
    do i=1, SYMBL
        do j=1, SYMBL
            X(i,j) = cmplx(1.0, 0.0, kind(0d0))
        end do
    end do

    V = PPL(H,HE,X,Nsybl,Npath,PPLloop)

	!average othogonality of eiven vector
    NAISEKI(1,1) = cmplx(0.0,0.0,kind(0d0))
    do i=1, SYMBL-1
        do j=i+1, SYMBL
            !固有ベクトル群を１列のベクトルに格納
            do k=1, SYMBL
                U(k,1) = X(k,i)
            end do
            !内積を取る固有ベクトルを格納
            do k=1, SYMBL
                N(k,1) = X(k,j)
            end do

            !随伴行列
            call CAdjoint(U,UH,SYMBL,1)

            !内積の計算
            call CMultiply(UH,N,UHN,1,SYMBL,SYMBL,1)

            !計算した内積を足し合わせる
            call CAdd(NAISEKI,UHN,NAISEKI,1,1,1,1)
        end do
    end do

    call CAbs(NAISEKI,NAISEKI_TMP,1,1)
    sC2 = SYMBL*(SYMBL-1.0)/2.0
    AVGOTH = NAISEKI_TMP / sC2

    do j=1, Nsybl
        do i=1, Nsybl
            print *, AVGOTH
        end do
    end do

end program