module CALmod
    implicit none
contains
	!行列の掛け算を行う

    subroutine CMultiply(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
        integer A_ROW,A_COL,B_ROW,B_COL,i,j,k
        complex(kind(0d0)) A(:,:), B(:,:),C(:,:)

        if(A_COL.ne.B_ROW) then
            print *, "can't calculate (Multiply)"
            stop
        end if

        do i=1, A_ROW
            do j=1, B_COL
                C(i,j) = cmplx(0.0,0.0,kind(0d0))
                do k=1, A_COL
                    C(i,j) = C(i,j)+A(i,k)*B(k,j)
                end do
            end do
        end do

    end subroutine

    !行列の和を取る

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

    !行列の差を取る

    subroutine CSubtract(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
        integer A_ROW,A_COL,B_ROW,B_COL,i,j
        complex(kind(0d0)) A(:,:),B(:,:),C(:,:)

        if((A_ROW.ne.B_ROW).or.(A_COL.ne.B_COL)) then
            print *, "can't calculate (Subtract)"
            stop
        end if

        do i=1, A_ROW
            do j=1, A_COL
                C(i,j) = A(i,j) - B(i,j)
            end do
        end do 

    end subroutine

    !随伴行列を取る

    subroutine CAdjoint(A,AH,A_ROW,A_COL)
        integer A_ROW,A_COL,i,j
        complex(kind(0d0)) A(:,:),AH(:,:)

        do i=1, A_ROW
            do j=1, A_COL
                AH(j,i) = conjg(A(i,j))
            end do
        end do

    end subroutine

    !処理Iを加える

    subroutine ProcI(A,AI,A_ROW,A_COL)
        integer A_ROW,A_COL,i,j
        complex(kind(0d0)) A(:,:),AI(:,:)

        do i=1, A_COL
            do j=1, A_ROW
                AI(A_ROW-j+1,i) = conjg(A(j,i))
            end do
        end do

    end subroutine

    !処理Jを加える

    subroutine ProcJ(A,AJ,A_ROW,A_COL,PATH)
        integer A_ROW,A_COL,PATH,i,j
        complex(kind(0d0)) A(:,:),AJ(:,:)

        do i=1, A_COL
            do j=PATH, A_ROW-(PATH-1)
                AJ(A_ROW-(PATH-1)-j+1,i) = conjg(A(j,i))
            end do
        end do

    end subroutine

    !配列を代入する

    subroutine CSubstitute(A,B,A_ROW,A_COL)
        integer A_ROW,A_COL,i,j
        complex(kind(0d0)) A(:,:),B(:,:)

        do i=1, A_ROW
            do j=1, A_COL
                A(i,j) = B(i,j)
            end do
        end do

    end subroutine

    !正規化する

    subroutine CNormalize(A,A_ROW,A_COL)
        integer A_ROW,A_COL,i
        complex(kind(0d0)) A(:,:)
        double precision TMP

        if(A_COL.ne.1) then
            print *, "The number of column isn't one."
            stop
        end if

        !実部と虚部の二乗の和を計算
        TMP = 0.0
        do i=1, A_ROW
            TMP = TMP + real(A(i,1))**2 + aimag(A(i,1))**2
        end do

        TMP = sqrt(TMP)

        !各成分をTMPで割る
        do i=1, A_ROW
            A(i,1) = A(i,1) / TMP
        end do

    end subroutine

    !複素数の絶対値をとる

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

    !正規乱数を配列に設定する

    subroutine Setnoise(A,A_ROW,A_COL)
        integer A_ROW,A_COL,i,j
        complex(kind(0d0)) A(:,:),CNV

        do i=1, A_ROW
            CNV=cmplx(Normal(),Normal(),kind(0d0))

            do j=1, A_COL
                A(i,j) = CNV
            end do
        end do

    end subroutine

    !正規乱数を返す
    function normal()
        double precision :: m=0.0
        double precision :: var=1.0
        double precision :: x=0.0
        double precision :: y=0.0
        double precision :: pi=dacos(dble(-1))
        double precision :: r1=0.0
        double precision :: r2=0.0
        double precision normal
        double precision :: r=0.0
        integer,save :: cnt=0

        integer seedsize,i
        integer,allocatable :: seed1(:)
        integer,allocatable :: seed2(:)

        !seedサイズの取得・動的配列の再宣言
        call random_seed(size=seedsize)
        allocate(seed1(seedsize))
        allocate(seed2(seedsize))

        if(cnt==0) then
            !初期シードを設定（独立になるようseed配列を書き換える）
            do i=1, seedsize
                seed1(i) = 1
            end do
            do i=1, seedsize
                seed2(i) = 1000000
            end do
            cnt = cnt + 1
        end if

        call random_seed(put=seed1)
        call random_number(r1)
        call random_seed(get=seed1)

        call random_seed(put=seed2)
        call random_number(r2)
        call random_seed(get=seed2)


        x = sqrt(-2.0*dlog(r1)*var)*dcos(2.0*pi*r2)+m
        y = sqrt(-2.0*dlog(r1)*var)*dsin(2.0*pi*r2)+m
        r = sqrt(x**2 + y**2)

        normal = x
    end function
end module CALmod