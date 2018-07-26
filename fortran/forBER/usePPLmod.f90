program usePPLmod
    use PPLmod
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

contains
    subroutine CAdjoint(A,AH,A_ROW,A_COL)
        integer A_ROW,A_COL,i,j
        complex(kind(0d0)) A(:,:),AH(:,:)

        do i=1, A_ROW
            do j=1, A_COL
                AH(j,i) = conjg(A(i,j))
            end do
        end do

    end subroutine

    subroutine CMultiply(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
        integer A_ROW,A_COL,B_ROW,B_COL,i,j,k
        complex(kind(0d0)) A(:,:), B(:,:),C(:,:)

        if(A_COL.ne.B_ROW) then
            print *, "can't calculate (Multiply)"
            stop
        end if

        do i=1, A_ROW
            do j=1, B_COL
                C(i,j) = cmplx(0.0,0.0, kind(0d0))
                do k=1, A_COL
                    C(i,j) = C(i,j)+A(i,k)*B(k,j)
                end do
            end do
        end do

    end subroutine

    subroutine CAdd(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
        integer A_ROW,A_COL,B_ROW,B_COL,i,j
        complex(kind(0d0)) A(:,:),B(:,:),C(:,:)

        if((A_ROW.ne.B_ROW).or.(A_COL.ne.B_COL)) then
            print *, "can't calculate (Add)"
            stop
        end if

        do i=1, A_ROW
            do j=1, A_COL
                C(i,j) = A(i,j) + B(i,j)
            end do
        end do 

    end subroutine

    subroutine CAbs(A,TMP,A_ROW,A_COL)
        integer i,A_ROW,A_COL
        complex(kind(0d0)) A(:,:)
        double precision TMP

        TMP = 0.0
        do i=1, A_ROW
            !実部と虚部の二乗の和を計算
            TMP = TMP + real(A(i,A_COL))**2 + aimag(A(i,A_COL))**2
        end do
        TMP = sqrt(TMP)
    end subroutine

end program