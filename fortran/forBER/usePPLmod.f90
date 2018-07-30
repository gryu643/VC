program usePPLmod
    use PPLmod
    use CALmod
    implicit none

    !declaration
    integer,parameter :: Nsybl=16
    integer,parameter :: Npath=8
    integer,parameter :: PPLloop=10000

    integer i,j,k
    complex(kind(0d0)) X(Nsybl,Nsybl)
    complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
    complex(kind(0d0)) V(Nsybl,Nsybl)
    complex(kind(0d0)) NAISEKI(1,1)
    complex(kind(0d0)) U(Nsybl,1)
    complex(kind(0d0)) N(Nsybl,1)
    complex(kind(0d0)) UH(1,Nsybl)
    complex(kind(0d0)) UHN(1,1)
    double precision NAISEKI_TMP
    double precision sC2
    double precision AVGOTH

    !initialize
    X(:,:)=(0.0,0.0)
    H(:,:)=(0.0,0.0)
    HE(:,:)=(0.0,0.0)
    V(:,:)=(0.0,0.0)
    NAISEKI(:,:)=(0.0,0.0)
    U(:,:)=(0.0,0.0)
    N(:,:)=(0.0,0.0)
    UH(:,:)=(0.0,0.0)
    UHN(:,:)=(0.0,0.0)
    NAISEKI_TMP=0.0
    sC2=0.0
    AVGOTH=0.0

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

    V = PPL(H,HE,X,Nsybl,Npath,PPLloop)

	!average othogonality of eiven vector
    NAISEKI(1,1) = cmplx(0.0,0.0,kind(0d0))
    do i=1, Nsybl-1
        do j=i+1, Nsybl
            !固有ベクトル群を１列のベクトルに格納
            do k=1, Nsybl
                U(k,1) = X(k,i)
            end do
            !内積を取る固有ベクトルを格納
            do k=1, Nsybl
                N(k,1) = X(k,j)
            end do

            !随伴行列
            call CAdjoint(U,UH,Nsybl,1)

            !内積の計算
            call CMultiply(UH,N,UHN,1,Nsybl,Nsybl,1)

            !計算した内積を足し合わせる
            call CAdd(NAISEKI,UHN,NAISEKI,1,1,1,1)
        end do
    end do

    call CAbs(NAISEKI,NAISEKI_TMP,1,1)
    sC2 = Nsybl*(Nsybl-1.0)/2.0
    AVGOTH = NAISEKI_TMP / sC2

    print *, AVGOTH


end program