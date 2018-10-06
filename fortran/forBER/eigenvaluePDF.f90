program eigenvaluePDF
    use CALmod
    use PPLmod
    implicit none

    !ifdef
    logical,parameter :: APPLY_PPL=.FALSE.

    !run time declaration
    integer t1, t2, t_rate, t_max, diff

    !declatation
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=2
    integer,parameter :: PPLloop=500
    integer,parameter :: trial=100000

    double precision Ampd(Npath,1)
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
    complex(kind(0d0)) Xppl(Nsybl,Nsybl)
    complex(kind(0d0)) HH(Nsybl,Nsybl+Npath-1)
    complex(kind(0d0)) HHH(Nsybl,Nsybl)
    complex(kind(0d0)) V(Nsybl,Nsybl)
    double precision Eig(1,Nsybl)
    complex(kind(0d0)) Eig_diag(1,Nsybl)

    double precision,parameter :: stride=0.01d0
    double precision output
    double precision,allocatable :: result(:,:)
    double precision,allocatable :: lambda(:)
    integer i,j,loop
    double precision start,end
    integer rank
    double precision avg
    double precision avgdB
    double precision x1
    double precision x2
    double precision lambdaP(1,Nsybl+1)

    !initialization
    Ampd(:,:)=0.0d0
    Cpath(:,:)=(0.0d0,0.0d0)
    H(:,:)=(0.0d0,0.0d0)
    HE(:,:)=(0.0d0,0.0d0)
    Xppl(:,:)=(0.0d0,0.0d0)
    HH(:,:)=(0.0d0,0.0d0)
    HHH(:,:)=(0.0d0,0.0d0)
    V(:,:)=(0.0d0,0.0d0)
    Eig(:,:)=0.0d0
    Eig_diag(:,:)=(0.0d0,0.0d0)
    avg=0.0d0
    lambdaP=0.0d0

    output=0.0d0
    start=0.0d0
    end=15.0d0
    !ファイル出力の要素数
    rank=nint((end-start)/stride)+1

    !allocate
    allocate(result(Nsybl+1,rank))
    allocate(lambda(rank))

    !allocate initialize
    result=0.0d0
    lambda=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open(7,file='evPDFdecomp_s32p2L1_0.0001.csv', status='replace')
    open(8,file='evPDFdecomp_s32p2L2_0.0001.csv', status='replace')
    open(9,file='evPDFdecomp_s32p2L3_0.0001.csv', status='replace')
    open(10,file='evPDFdecomp_s32p2L4_0.0001.csv', status='replace')
    open(11,file='evPDFdecomp_s32p2L5_0.0001.csv', status='replace')
    open(12,file='evPDFdecomp_s32p2L6_0.0001.csv', status='replace')
    open(13,file='evPDFdecomp_s32p2L7_0.0001.csv', status='replace')
    open(14,file='evPDFdecomp_s32p2L8_0.0001.csv', status='replace')
    open(15,file='evPDFdecomp_s32p2L9_0.0001.csv', status='replace')
    open(16,file='evPDFdecomp_s32p2L10_0.0001.csv', status='replace')
    open(17,file='evPDFdecomp_s32p2L11_0.0001.csv', status='replace')
    open(18,file='evPDFdecomp_s32p2L12_0.0001.csv', status='replace')
    open(19,file='evPDFdecomp_s32p2L13_0.0001.csv', status='replace')
    open(20,file='evPDFdecomp_s32p2L14_0.0001.csv', status='replace')
    open(21,file='evPDFdecomp_s32p2L15_0.0001.csv', status='replace')
    open(22,file='evPDFdecomp_s32p2L16_0.0001.csv', status='replace')
    open(23,file='evPDFdecomp_s32p2L17_0.0001.csv', status='replace')
    open(24,file='evPDFdecomp_s32p2L18_0.0001.csv', status='replace')
    open(25,file='evPDFdecomp_s32p2L19_0.0001.csv', status='replace')
    open(26,file='evPDFdecomp_s32p2L20_0.0001.csv', status='replace')
    open(27,file='evPDFdecomp_s32p2L21_0.0001.csv', status='replace')
    open(28,file='evPDFdecomp_s32p2L22_0.0001.csv', status='replace')
    open(29,file='evPDFdecomp_s32p2L23_0.0001.csv', status='replace')
    open(30,file='evPDFdecomp_s32p2L24_0.0001.csv', status='replace')
    open(31,file='evPDFdecomp_s32p2L25_0.0001.csv', status='replace')
    open(32,file='evPDFdecomp_s32p2L26_0.0001.csv', status='replace')
    open(33,file='evPDFdecomp_s32p2L27_0.0001.csv', status='replace')
    open(34,file='evPDFdecomp_s32p2L28_0.0001.csv', status='replace')
    open(35,file='evPDFdecomp_s32p2L29_0.0001.csv', status='replace')
    open(36,file='evPDFdecomp_s32p2L30_0.0001.csv', status='replace')
    open(37,file='evPDFdecomp_s32p2L31_0.0001.csv', status='replace')
    open(38,file='evPDFdecomp_s32p2L32_0.0001.csv', status='replace')
    open(39,file='evPDFdecomp_s32p2_0.0001.csv', status='replace')

    !implementation
    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do


	do loop=1, trial
        do i=1, Npath
            Cpath(i,1) = cmplx(normal(),normal(),kind(0d0))/sqrt(2.0d0)*Ampd(i,1)
        end do

        !set H
        do i=1, Nsybl
            do j=1, Npath
                H(i+j-1,i) = Cpath(j,1)
            end do
        end do

        !set HE
        do i=1, Nsybl+Npath-1
            do j=1, Npath
                HE(i+j-1,i) = Cpath(j,1)
            end do
        end do

        !set HH
        call CAdjoint(H,HH,Nsybl+Npath-1,Nsybl)

        !set HHH
        call CMultiply(HH,H,HHH,Nsybl,Nsybl+Npath-1,Nsybl+Npath-1,Nsybl)

        if(APPLY_PPL) then
            !set Xppl
            do i=1, Nsybl
                do j=1, Nsybl
                    Xppl(i,j) = cmplx(1.0d0, 0.0d0, kind(0d0))
                end do
            end do

            call PPL(H,HE,Xppl,Eig,Nsybl,Npath,PPLloop)
        else
            call CSubstitute(V,HHH,Nsybl,Nsybl)
!            call decomp_zheev(Nsybl,V,Eig)
!            call decomp_zheevd(Nsybl,V,Eig)
!            call decomp_zgeev(Nsybl,V,Eig)
            call decomp_zhpev(Nsybl,V,Eig)
        endif

        call sort(Eig,Nsybl)

		do i=1, Nsybl
            output = Eig(1,i)
            do j=1, rank
                x1 = start+stride*(j-1)
                x2 = start+stride*j
                if((output.ge.x1).and.(output.lt.x2)) then
                    result(i,j) = result(i,j) + 1.0d0
                    exit
                endif
            end do
        end do
	end do
    
    !calculate lambda 1~Nsybl
    do i=1, rank
        lambda(i) = start + stride*dble(i-1)
    end do

    !calculate p(lambda)
    do i=1, Nsybl
        do j=1, rank
            result(i,j) = result(i,j)/trial
        end do
    end do
    
    !calculate p(lambda) average lambda
    do i=1, rank
        do j=1, Nsybl
            result(Nsybl+1,i) = result(Nsybl+1,i) + result(j,i)
        end do
    end do
    result(Nsybl+1,:) = result(Nsybl+1,:)/(trial*dble(Nsybl))

    !normalize p(lambda)
    do i=1, Nsybl+1
        do j=1, rank
            lambdaP(1,i) = lambdaP(1,i) + lambda(j)*result(i,j)
        end do
    end do
    
    do i=1, Nsybl+1
        do j=1, rank
            result(i,j) = result(i,j)/lambdaP(1,i)
        end do
    end do

    do j=1, rank
        write(7,*) lambda(j), ',', result(1,j)
        write(8,*) lambda(j), ',', result(2,j)
        write(9,*) lambda(j), ',', result(3,j)
        write(10,*) lambda(j), ',', result(4,j)
        write(11,*) lambda(j), ',', result(5,j)
        write(12,*) lambda(j), ',', result(6,j)
        write(13,*) lambda(j), ',', result(7,j)
        write(14,*) lambda(j), ',', result(8,j)
        write(15,*) lambda(j), ',', result(9,j)
        write(16,*) lambda(j), ',', result(10,j)
        write(17,*) lambda(j), ',', result(11,j)
        write(18,*) lambda(j), ',', result(12,j)
        write(19,*) lambda(j), ',', result(13,j)
        write(20,*) lambda(j), ',', result(14,j)
        write(21,*) lambda(j), ',', result(15,j)
        write(22,*) lambda(j), ',', result(16,j)
        write(23,*) lambda(j), ',', result(17,j)
        write(24,*) lambda(j), ',', result(18,j)
        write(25,*) lambda(j), ',', result(19,j)
        write(26,*) lambda(j), ',', result(20,j)
        write(27,*) lambda(j), ',', result(21,j)
        write(28,*) lambda(j), ',', result(22,j)
        write(29,*) lambda(j), ',', result(23,j)
        write(30,*) lambda(j), ',', result(24,j)
        write(31,*) lambda(j), ',', result(25,j)
        write(32,*) lambda(j), ',', result(26,j)
        write(33,*) lambda(j), ',', result(27,j)
        write(34,*) lambda(j), ',', result(28,j)
        write(35,*) lambda(j), ',', result(29,j)
        write(36,*) lambda(j), ',', result(30,j)
        write(37,*) lambda(j), ',', result(31,j)
        write(38,*) lambda(j), ',', result(32,j)
        write(39,*) lambda(j), ',', result(33,j)
    end do

    !file close
    do i=1+6, Nsybl+7
        close(i)
    end do

    !time measurement end
    call system_clock(t2,t_rate,t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program