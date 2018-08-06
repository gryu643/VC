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

    !normal random number
    function normal()
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
        y = sqrt(-2.0d0*dlog(r1)*var)*dsin(2.0d0*pi*r2)+m
        r = sqrt(x**2 + y**2)

        normal = x
    end function

    !random number
	function rand()
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
    
    !return eigenvalue decomposition
    subroutine zdiag(N,A,ev)
        !sikinote
        implicit none
        integer::N
        double precision::ev(1:N)
        complex(kind(0d0))::A(1:N,1:N)
       
        integer::lda,lwork,liwork,lrwork,info
        double precision,allocatable::rwork(:)
        complex(kind(0d0)),allocatable::work(:)
        complex(kind(0d0))::qw(1:3)
        double precision::qrw(1:3)
        integer,allocatable::iwork(:)
        integer::qiw(1:3)
        character(1)::job,uplo
       
        job="V"
        uplo="U"
        lda=N
       
        call zheevd(job,uplo,N,0,lda,0,qw,-1,qrw,-1,qiw,-1,info)
        if(info.ne.0)then
           write(6,'(A)')"    program stop @zheevd"
           write(6,'(A,i0)')"    info --> ",info
           stop
        endif
       
        lrwork=idint(qrw(1))+1
        lwork=idint(dble(qw(1)))+1
        liwork=qiw(1)+1
        allocate(work(1:lwork),iwork(1:liwork),rwork(1:lrwork))
        work(1:lwork)=0.d0
        iwork(1:liwork)=0
        rwork(1:lrwork)=0.d0
       
        !diagonalise.
        call zheevd(job,uplo,N,A,lda,ev,work,lwork,rwork,lrwork,iwork,liwork,info)
        if(info.ne.0)then
           write(6,'(A)')"    program stop @zheevd"
           write(6,'(A,i0)')"    info --> ",info
           stop
        endif
       
        deallocate(work,iwork)
       
        return
      end subroutine zdiag

    subroutine diag(N,A,Ev)
        ! sikinote
        !date      : 2015/07/07
        !            2015/08/21   
        !developer : sikino & fernandeskun
        implicit none
        integer,intent(in)::N
        complex(kind(0d0)),intent(inout)::A(1:N,1:N)
        complex(kind(0d0)),intent(out)::Ev(1:N)

        integer::ilo,ihi,info,lwork,turn(1:N),tmp,i
        double precision::scale(1:N),rwork(1:N)
        complex(kind(0d0))::tau(1:N-1),w(1:N),z(1:N,1:N),Q(1:N,1:N),vr(1:N,1:N),tw(1:3)
        complex(kind(0d0)),allocatable::work(:)
        
        tau(1:N-1)=dcmplx(0d0,0d0)
        w(1:N)=dcmplx(0d0,0d0)
        z(1:N,1:N)=dcmplx(0d0,0d0)
        Q(1:N,1:N)=dcmplx(0d0,0d0)
        vr(1:N,1:N)=dcmplx(0d0,0d0)
        tw(1:3)=dcmplx(0d0,0d0)
        Ev(1:N)=dcmplx(0d0,0d0)
        
        !Equilibrate matrix A to equilibrated matrix A' to improve accuracy.  
        !            i   i io  i   o    o     o     o 
        call zgebal('P', N, A, N, ilo, ihi, scale, info)
        if(info.ne.0)then
            write(6,'(A,i0)')" At zgebal error, info --> ",info
            write(6,'(A)')" Program stop"
            stop
        endif
        
        !Size Query
        call zgehrd(N, ilo, ihi, A, N, tau, tw, -1, info)
        lwork=nint(dble(tw(1)))
        allocate(work(1:lwork)); work=dcmplx(0d0,0d0)
        
        !Degenerate matrix A to upper Hessenberg matrix H.   
        !           i   i    i  io  i   o    i      i      o
        call zgehrd(N, ilo, ihi, A, N, tau, work, lwork, info)
        if(info.ne.0)then
            write(6,'(A,i0)')" At zgehrd error, info --> ",info
            write(6,'(A)')" Program stop"
            stop
        endif
        deallocate(work)

        Q=a
        !Size Query
        call zunghr(N, ilo, ihi, Q, N, tau, tw, -1, info)
        lwork=nint(dble(tw(1)))
        allocate(work(1:lwork)); work=dcmplx(0d0,0d0)

        !Make complex unitary matrix Q from upper Hessenberg matrix H.
        !           i   i    i  io  i   i    i      i      o
        call zunghr(N, ilo, ihi, Q, N, tau, work, lwork, info)
        if(info.ne.0)then
            write(6,'(A,i0)')" At zunghr error, info --> ",info
            write(6,'(A)')" Program stop"
            stop
        endif
        deallocate(work)
        
        z=Q
        !Size Query
        call zhseqr('S', 'V', N, ilo, ihi, A, N, Ev, z, N, tw, -1, info)
        lwork=nint(dble(tw(1)))
        allocate(work(1:lwork)); work=dcmplx(0d0,0d0)

        !Get eigenvalue of upper Hessenberg matrix H and Get Schur vector.
        !                     i   i    i  io  i   o  o  i   i      i      o  
        call zhseqr('S', 'V', N, ilo, ihi, A, N, Ev, z, N, work, lwork, info)
        if(info.ne.0)then
            write(6,'(A,i0)')" At zhseqr error, info --> ",info
            write(6,'(A)')" Program stop"
            stop
        endif
        deallocate(work)

        !Get right eigenvector X from upper triangular matrix T. 
        allocate(work(1:2*N))  
        vr=z
        !                        i  i  i         o  i  i   o   i      i      i
        call ztrevc('R', 'B', 0, N, A, N, 0, 1, vr, N, N, tmp, work, rwork, info)
        if(info.ne.0)then
            write(6,'(A,i0)')" At zhseqr error, info --> ",info
            write(6,'(A)')" Program stop"
            stop
        endif
        deallocate(work)
        
        !Transrate right eigenvector X of Equilibrated matrix A' to right eigenvector of matrix A
        !                     i   i    i     i    i   o  i   o 
        call zgebak('P', 'R', N, ilo, ihi, scale, N, vr, N, info)
        if(info.ne.0)then
            write(6,'(A,i0)')" At zhseqr error, info --> ",info
            write(6,'(A)')" Program stop"
            stop
        endif
            
        A=vr

        !swap Eigenvectol as same arrangement for Eigenvalue
        call sortdp2(N,Ev,turn)

        Q=A
        do i=1,N
            tmp=turn(i)
            A(1:N,i)=Q(1:N,tmp)
        enddo
        return

    !sort Eigenvalue of real part from small to big.
    contains
        subroutine sortdp2(N,data,turn)
            implicit none
            integer::i,ti,j,N,turn(1:N)
            complex(kind(0d0))::data(1:N),tmp

            do i=1,N
            turn(i)=i
            enddo

            do i=1,N-1
            do j=i+1,N
                if(dble(data(i)) > dble(data(j)))then
                    tmp=data(i)
                    data(i)=data(j)
                    data(j)=tmp

                    ti=turn(i)
                    turn(i)=turn(j)
                    turn(j)=ti
                end if
            end do
            end do

            return
        end subroutine sortdp2
    end subroutine diag

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

        !-- return
        do i=1, Nsybl
            Eig(1,i) = real(w(i))
        end do

    end subroutine decomp_zgeev
    
end module CALmod