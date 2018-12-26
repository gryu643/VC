module CALmod
    implicit none

    double precision :: pi=dacos(dble(-1))
contains
    subroutine CMultiply(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
        !--------------------------------------------------------------------------!
        !complex multiplication
        !--------------------------------------------------------------------------!

        integer A_ROW,A_COL,B_ROW,B_COL,i,j,k
        complex(kind(0d0)) A(:,:), B(:,:),C(:,:)
        complex(kind(0d0)) db
        C=(0.0d0,0.0d0)

        do j=1, B_COL
            do k=1, B_ROW
                db = B(k,j)
                do i=1, A_ROW
                    C(i,j) = C(i,j)+A(i,k)*db
                end do
            end do
        end do

    end subroutine

    subroutine RMultiply(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
        !--------------------------------------------------------------------------!
        !real number multiplication
        !--------------------------------------------------------------------------!

        integer A_ROW,A_COL,B_ROW,B_COL,i,j,k
        double precision A(:,:), B(:,:),C(:,:)
        double precision db

        do j=1, B_COL
            do k=1, B_ROW
                db = B(k,j)
                do i=1, A_ROW
                    C(i,j) = C(i,j)+A(i,k)*db
                end do
            end do
        end do

    end subroutine

    subroutine CAdd(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
        !--------------------------------------------------------------------------!
        !complex Addition
        !--------------------------------------------------------------------------!

        integer A_ROW,A_COL,B_ROW,B_COL,i,j
        complex(kind(0d0)) A(:,:),B(:,:),C(:,:)

        do j=1, A_COL
            do i=1, A_ROW
                C(i,j) = A(i,j) + B(i,j)
            end do
        end do 

    end subroutine

    subroutine CSubtract(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
        !--------------------------------------------------------------------------!
        !complex subtraction
        !--------------------------------------------------------------------------!

        integer A_ROW,A_COL,B_ROW,B_COL,i,j
        complex(kind(0d0)) A(:,:),B(:,:),C(:,:)

        do j=1, A_COL
            do i=1, A_ROW
                C(i,j) = A(i,j) - B(i,j)
            end do
        end do 

    end subroutine

    subroutine CAdjoint(A,AH,A_ROW,A_COL)
        !--------------------------------------------------------------------------!
        !calculate Hermitian conjugate
        !--------------------------------------------------------------------------!

        integer A_ROW,A_COL,i,j
        complex(kind(0d0)) A(:,:),AH(:,:)

        do j=1, A_COL
            do i=1, A_ROW
                AH(j,i) = conjg(A(i,j))
            end do
        end do

    end subroutine

    subroutine ProcI(A,AI,A_ROW,A_COL)
        !--------------------------------------------------------------------------!
        !do processing I(subroutine used in PPL)
        !--------------------------------------------------------------------------!

        integer A_ROW,A_COL,i,j
        complex(kind(0d0)) A(:,:),AI(:,:)

        do i=1, A_COL
            do j=1, A_ROW
                AI(A_ROW-j+1,i) = conjg(A(j,i))
            end do
        end do

    end subroutine

    subroutine ProcJ(A,AJ,A_ROW,A_COL,PATH)
        !--------------------------------------------------------------------------!
        !do processing J(subroutine used in PPL)
        !--------------------------------------------------------------------------!

        integer A_ROW,A_COL,PATH,i,j
        complex(kind(0d0)) A(:,:),AJ(:,:)

        do i=1, A_COL
            do j=PATH, A_ROW-(PATH-1)
                AJ(A_ROW-(PATH-1)-j+1,i) = conjg(A(j,i))
            end do
        end do

    end subroutine

    subroutine CSubstitute(A,B,A_ROW,A_COL)
        !--------------------------------------------------------------------------!
        !substitute complex
        !--------------------------------------------------------------------------!

        integer A_ROW,A_COL,i,j
        complex(kind(0d0)) A(:,:),B(:,:)

        do j=1, A_COL
            do i=1, A_ROW
                A(i,j) = B(i,j)
            end do
        end do

    end subroutine

    subroutine CNormalize(A,A_ROW,A_COL)
        !--------------------------------------------------------------------------!
        !Normalize complex
        !--------------------------------------------------------------------------!

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

    subroutine CAbs(A,TMP,A_ROW,A_COL)
        !--------------------------------------------------------------------------!
        !return the absolute value
        !--------------------------------------------------------------------------!

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

    function normal()
        !------------------------------------------------------------------------!
        !normal random number(mean=0,variance=1)
        !------------------------------------------------------------------------!

        double precision :: m=0.0d0
        double precision :: var=1.0d0
        double precision :: x=0.0d0
        double precision :: y=0.0d0
        double precision :: pi=dacos(dble(-1))
        double precision :: r1=0.0d0
        double precision :: r2=0.0d0
        double precision normal
        double precision :: r=0.0d0
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
                call system_clock(count=seed1(i))
                seed1(i) = 1
            end do
            do i=1, seedsize
            call system_clock(count=seed2(i))
                seed2(i) = seed2(i) - 100000
            end do
            cnt = cnt + 1
        end if

        call random_seed(put=seed1)
        call random_number(r1)
        call random_seed(get=seed1)

        call random_seed(put=seed2)
        call random_number(r2)
        call random_seed(get=seed2)


        x = sqrt(-2.0d0*dlog(r1)*var)*dcos(2.0d0*pi*r2)+m
!        y = sqrt(-2.0d0*dlog(r1)*var)*dsin(2.0d0*pi*r2)+m
!        r = sqrt(x**2 + y**2)

        normal = x
    end function

	function rand()
        !----------------------------------------------------------------------------!
        !uniform random number[0,1)
        !----------------------------------------------------------------------------!

		!declaration
		integer seedsize,i
		integer,allocatable :: seed(:)
		integer,save :: rand_cnt=0
		double precision rand

		!initialization
		rand=0.0

		call random_seed(size=seedsize)
		allocate(seed(seedsize))

		if(rand_cnt==0) then
			do i=1, seedsize
                call system_clock(count=seed(i))
			end do
		end if

		call random_number(rand)
    end function

    function orthogonal(X,Nsybl)
        !--------------------------------------------------------------------------!
        !return orthogonality
        !--------------------------------------------------------------------------!

        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) X(Nsybl,Nsybl)

        !-- declaration
        integer i,j,k
        complex(kind(0d0)) NAISEKI(1,1)
        complex(kind(0d0)) V1(Nsybl)
        complex(kind(0d0)) V2(Nsybl)
        double precision NAISEKI_TMP
        double precision orthogonal

        !-- initialization
        NAISEKI=(0.0d0,0.0d0)
        V1=(0.0d0,0.0d0)
        V2=(0.0d0,0.0d0)
        NAISEKI_TMP=0.0d0
        orthogonal=0.0d0

        !-- implementation
			!average othogonality of eiven vector
			NAISEKI(1,1) = cmplx(0.0,0.0,kind(0d0))
			do i=1, Nsybl-1
				!固有ベクトル群を１列のベクトルに格納
				do k=1, Nsybl
					V1(k) = X(k,i)
				end do
				do j=i+1, Nsybl
					!内積を取る固有ベクトルを格納
					do k=1, Nsybl
						V2(k) = X(k,j)
					end do
					NAISEKI(1,1) = NAISEKI(1,1) + abs(dot_product(V1,V2))
				end do
			end do

			call CAbs(NAISEKI,NAISEKI_TMP,1,1)
            orthogonal = NAISEKI_TMP
    end function orthogonal

    subroutine PPLEV(X,HHHX,LAMBDA,Nsybl)
        !----------------------------------------------------------------!
        !calculate eigenvalue(subroutine used in PPL)
        !----------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) X(Nsybl,Nsybl)
        complex(kind(0d0)) HHHX(Nsybl,Nsybl)
        complex(kind(0d0)) LAMBDA(Nsybl,1)

        !-- declaration
        integer i,j
        complex(kind(0d0)) D_X(Nsybl,1)
        complex(kind(0d0)) D_HHHX(Nsybl,1)
        double precision D_XAbs
        double precision D_HHHXAbs
        
        !-- initialization
        D_X=(0.0d0,0.0d0)
        D_HHHX=(0.0d0,0.0d0)
        D_XAbs=0.0d0
        D_HHHXAbs=0.0d0

        !-- implementation
        do i=1, Nsybl
            do j=1,Nsybl
                D_X(j,1) = X(j,i)
                D_HHHX(j,1) = HHHX(j,i)
            end do

            call CAbs(D_X,D_XAbs,Nsybl,1)
            call CAbs(D_HHHX,D_HHHXAbs,Nsybl,1)

            LAMBDA(i,1) = cmplx(D_HHHXAbs/D_XAbs,0.0d0, kind(0d0))
        end do
    end subroutine PPLEV

    subroutine PPLSubPart(X,LAMBDA,SUB_PART,Nsybl)
        !-----------------------------------------------------------------------!
        !calculate Subpart(subroutine used in PPL)
        !-----------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) X(Nsybl,Nsybl)
        complex(kind(0d0)) SUB_PART(Nsybl,Nsybl)
        complex(kind(0d0)) LAMBDA(Nsybl,1)

        !-- declaration
        integer i,k
		complex(kind(0d0)) Xn(Nsybl,1)
		complex(kind(0d0)) U(Nsybl,1)
		complex(kind(0d0)) UH(1,Nsybl)
		complex(kind(0d0)) LUUH(Nsybl,Nsybl)
		complex(kind(0d0)) LUUH_SET(Nsybl,Nsybl)
		complex(kind(0d0)) LU(Nsybl,1)
		complex(kind(0d0)) LUUHXn(Nsybl,1)

        !-- initialization
		Xn=(0.0,0.0)
		U=(0.0,0.0)
		UH=(0.0,0.0)
		LUUH=(0.0,0.0)
		LUUH_SET=(0.0,0.0)
		LU=(0.0,0.0)
		LUUHXn=(0.0,0.0)

        !-- implementation
        LUUH_SET=cmplx(0.0d0,0.0d0,kind(0d0))
        do i=2, Nsybl
            !収束する固有ベクトル(Nsybl,1)
            do k=1, Nsybl
                Xn(k,1) = X(k,i)
            end do

            !減算する固有ベクトル(Nsybl,1)
            do k=1, Nsybl
                U(k,1) = X(k,i-1)
            end do

            !固有ベクトルの随伴行列(1,Nsybl)
            call CAdjoint(U,UH,Nsybl,1)

            !λ*U(Nsybl,1)
            do k=1, Nsybl
                LU(k,1) = LAMBDA(i-1,1)*U(k,1)
            end do

            !LU*UH(Nsybl,Nsybl)
            call CMultiply(LU,UH,LUUH,Nsybl,1,1,Nsybl)

            !λUUHの集合を格納
            call CAdd(LUUH_SET,LUUH,LUUH_SET,Nsybl,Nsybl,Nsybl,Nsybl)

            !LUUH_SET*Xn(Nsybl,1) iの次ループで足し合わせる
            call CMultiply(LUUH_SET,Xn,LUUHXn,Nsybl,Nsybl,Nsybl,1)

            !減算部の格納
            do k=1, Nsybl
                SUB_PART(k,i) = LUUHXn(k,1)
            end do
        end do
    end subroutine PPLSubPart

    subroutine PPLNormalize(arSUB,X,Nsybl)
        !-----------------------------------------------------------!
        !Normalize Eigen vector(subroutine used in PPL)
        !-----------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) arSUB(Nsybl,Nsybl)
        complex(kind(0d0)) X(Nsybl,Nsybl)

        !-- declaration
        integer i,j
        complex(kind(0d0)) NORM(Nsybl,1)

        !-- initialization
        NORM=(0.0d0,0.0d0)

        !-- implementation
        do i=1, Nsybl
            !固有ベクトル群を1列のベクトルに格納
            do j=1, Nsybl
                NORM(j,1) = arSUB(j,i)
            end do

            !正規化
            call CNormalize(NORM,Nsybl,1)

            !正規化したベクトルをXに格納
            do j=1, Nsybl
                X(j,i) = NORM(j,1)
            end do
        end do
    end subroutine PPLNormalize

    subroutine decomp_zheevd(Nsybl,V,Eig)
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) V(Nsybl,Nsybl)
        double precision Eig(1,Nsybl)

        !-- declaration
        integer :: i,ifail,info,lda,liwork,lrwork,lwork,n
        character(1) :: job, uplo
        complex(kind(0d0)),allocatable :: work(:)
        double precision,allocatable :: rwork(:),w(:)
        integer,allocatable :: iwork(:)

        !-- initialization
        n = Nsybl
        lda = n
        liwork = 5*n + 3
        lrwork = 2*n*n + 5*n + 1
        lwork = n*(n+2)
        job = 'V'
        uplo = 'L'

        !-- allocate
        allocate(work(lwork),rwork(lrwork),w(n),iwork(liwork))

        !-- implementation
        ! calculate all the eigenvalues and eigenvectors
        call zheevd(job,uplo,n,V,lda,w,work,lwork,rwork,lrwork,iwork, &
        liwork,info)

        !-- return
        do i=1, Nsybl
            Eig(1,i) = w(i)
        end do

    end subroutine decomp_zheevd

    subroutine decomp_zheev(Nsybl,V,Eig)
        !--------------------------------------------------------------------------!
        !eigenvalue decomposition(using zheev nagfor library)
        !--------------------------------------------------------------------------!

        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) V(Nsybl,Nsybl)
        double precision Eig(1,Nsybl)

        !-- declaration
        integer :: i,ifail,info,lda,lrwork,lwork,n
        character(1) :: job, uplo
        complex(kind(0d0)),allocatable :: work(:)
        double precision,allocatable :: rwork(:),w(:)

        !-- initialization
        n = Nsybl
        lda = n
        lrwork = 3*n-2
        lwork = n*(n+2)
        job = 'V'
        uplo = 'L'

        !-- allocate
        allocate(work(lwork),rwork(lrwork),w(n))

        !-- implementation
        ! calculate all the eigenvalues and eigenvectors
        call zheev(job,uplo,n,V,lda,w,work,lwork,rwork,info)

        !-- return
        do i=1, Nsybl
            Eig(1,i) = w(i)
        end do

    end subroutine decomp_zheev

    subroutine decomp_zgeev(Nsybl,V,Eig)
        !--------------------------------------------------------------------------!
        !eigenvalue decomposition(using zgeev nagfor library)
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) V(Nsybl,Nsybl)
        double precision Eig(1,Nsybl)

        !-- declaration
        integer,parameter :: nb=64
        integer :: i,ifail,info,lda,ldvr,lwork,n
        complex(kind(0d0)) :: dummy(1,1)
        complex(kind(0d0)),allocatable :: work(:),vr(:,:),w(:)
        double precision,allocatable :: rwork(:)

        !-- initialization
        n = Nsybl
        lda = n
        ldvr = n

        !-- allocate
        allocate(rwork(2*n),w(n),vr(ldvr,n))

        !-- implementation
        ! use routine workspace query to get optimal workspace.
        lwork = -1
        call zgeev('No left vectors','Vectors(right)',n,V,lda,w,dummy,1, &
        vr,ldvr,dummy,lwork,rwork,info)

        !Make sure that there is enough workspace for block size nb.
        lwork = max((nb+1)*n, nint(real(dummy(1,1))))
        allocate(work(lwork))

        ! calculate all the eigenvalues and eigenvectors
        call zgeev('No left vectors','Vectors(right)',n,V,lda,w,dummy,1, &
        vr,ldvr,work,lwork,rwork,info)

        if(info==0) then
            !-- return
            do i=1, Nsybl
                Eig(1,i) = real(w(i))
            end do
        else
            print *, 'Failure in ZGEEV. info =', info
        end if

    end subroutine decomp_zgeev

    subroutine decomp_zhpev(Nsybl,V,Eig)
        !--------------------------------------------------------------------------!
        !eigenvalue decomposition(using zhpev nagfor library)
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        integer Nsybl
        complex(kind(0d0)) V(Nsybl,Nsybl)
        double precision Eig(1,Nsybl)

        !-- declaration
        integer :: i,j,info,n
        character(1),parameter :: uplo='U'
        complex(kind(0d0)) :: dummy(1,1)
        complex(kind(0d0)),allocatable :: work(:),ap(:)
        double precision,allocatable :: rwork(:),w(:)

        !-- initialization
        n = Nsybl

        !-- allocate
        allocate(ap((n*(n+1))/2),rwork(3*n-2),w(n),work(2*n-1))

        !-- implementation
        do i=1, n
            do j=i,n
                ap(i+(j*(j-1))/2) = V(i,j)
            end do
        end do

        ! calculate all the eigenvalues and eigenvectors
        call zhpev('V', uplo, n, ap, w, V, n, work, rwork, info)

        if(info==0) then
            !-- return
            do i=1, Nsybl
                Eig(1,i) = real(w(i))
            end do
        else
            print *, 'Failure in ZGEEV. info =', info
        end if

    end subroutine decomp_zhpev

    subroutine sort(A,num)
        !--------------------------------------------------------------------------!
        !descending sort
        !--------------------------------------------------------------------------!

        implicit none

        !-- declaration
        integer i,j,num
        double precision A(1,num)
        double precision TMP

        !-- initialization
        TMP=0.0d0

        !-- implementation
        do i=1, num
            do j=i+1, num
                if(A(1,i)<A(1,j)) then
                    TMP = A(1,i)
                    A(1,i) = A(1,j)
                    A(1,j) = TMP
                endif
            end do
        end do

    end subroutine sort

    subroutine InverseMat(A,n)
        !--------------------------------------------------------------------------!
        !calculate inverse matrix
        !--------------------------------------------------------------------------!
!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: nin = 5, nout = 6
!      .. Local Scalars ..
       INTEGER                         :: i, ifail, info, lda, lwork, n
!      .. Local Arrays ..
       double precision, ALLOCATABLE :: work(:)
       double precision A(n,n)
       INTEGER, ALLOCATABLE            :: ipiv(:)
!      .. Executable Statements ..

!      Skip heading in data file
       lda = n
       lwork = 64*n
       ALLOCATE (work(lwork),ipiv(n))

!      Factorize A

!      The NAG name equivalent of dgetrf is f07adf
       CALL dgetrf(n,n,A,lda,ipiv,info)

       FLUSH (nout)
       IF (info==0) THEN

!         Compute inverse of A

!         The NAG name equivalent of dgetri is f07ajf
          CALL dgetri(n,A,lda,ipiv,work,lwork,info)

!         Print inverse

!         ifail: behaviour on error exit
!                =0 for hard exit, =1 for quiet-soft, =-1 for noisy-soft
          ifail = 0
!          CALL x04caf('General',' ',n,n,A,lda,'Inverse',ifail)

       ELSE
          WRITE (nout,*) 'The factor U is singular'
       END IF

    end subroutine InverseMat

    function euler(theta)
        !--------------------------------------------------------------------------!
        !return euler
        !--------------------------------------------------------------------------!
        implicit none

        !-- argument
        double precision theta

        !-- decralation
        complex(kind(0d0)) euler

        !-- initialization
        euler=(0.0d0,0.0d0)

        !-- implementation
        euler = cmplx(dcos(theta),dsin(theta),kind(0d0))
    end function euler
end module CALmod
