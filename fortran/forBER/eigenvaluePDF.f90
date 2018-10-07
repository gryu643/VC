program eigenvaluePDF
    use CALmod
    implicit none

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
    double precision,parameter :: start=0.0d0
    double precision,parameter :: end=20.0d0
    integer rank
    double precision avg
    double precision avgdB
    double precision x1
    double precision x2
    double precision lambdaP(1,Nsybl+1)
    integer OutPutRow

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
    OutPutRow=0
    output=0.0d0
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
    open(7,file='ev(s32p2L1)0.01.csv', status='replace')
    open(8,file='ev(s32p2L2)0.01.csv', status='replace')
    open(9,file='ev(s32p2L3)0.01.csv', status='replace')
    open(10,file='ev(s32p2L4)0.01.csv', status='replace')
    open(11,file='ev(s32p2L5)0.01.csv', status='replace')
    open(12,file='ev(s32p2L6)0.01.csv', status='replace')
    open(13,file='ev(s32p2L7)0.01.csv', status='replace')
    open(14,file='ev(s32p2L8)0.01.csv', status='replace')
    open(15,file='ev(s32p2L9)0.01.csv', status='replace')
    open(16,file='ev(s32p2L10)0.01.csv', status='replace')
    open(17,file='ev(s32p2L11)0.01.csv', status='replace')
    open(18,file='ev(s32p2L12)0.01.csv', status='replace')
    open(19,file='ev(s32p2L13)0.01.csv', status='replace')
    open(20,file='ev(s32p2L14)0.01.csv', status='replace')
    open(21,file='ev(s32p2L15)0.01.csv', status='replace')
    open(22,file='ev(s32p2L16)0.01.csv', status='replace')
    open(23,file='ev(s32p2L17)0.01.csv', status='replace')
    open(24,file='ev(s32p2L18)0.01.csv', status='replace')
    open(25,file='ev(s32p2L19)0.01.csv', status='replace')
    open(26,file='ev(s32p2L20)0.01.csv', status='replace')
    open(27,file='ev(s32p2L21)0.01.csv', status='replace')
    open(28,file='ev(s32p2L22)0.01.csv', status='replace')
    open(29,file='ev(s32p2L23)0.01.csv', status='replace')
    open(30,file='ev(s32p2L24)0.01.csv', status='replace')
    open(31,file='ev(s32p2L25)0.01.csv', status='replace')
    open(32,file='ev(s32p2L26)0.01.csv', status='replace')
    open(33,file='ev(s32p2L27)0.01.csv', status='replace')
    open(34,file='ev(s32p2L28)0.01.csv', status='replace')
    open(35,file='ev(s32p2L29)0.01.csv', status='replace')
    open(36,file='ev(s32p2L30)0.01.csv', status='replace')
    open(37,file='ev(s32p2L31)0.01.csv', status='replace')
    open(38,file='ev(s32p2L32)0.01.csv', status='replace')
    open(39,file='ev(s32p2)0.01.csv', status='replace')

    !implementation
    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do


	do loop=1, trial
        do i=1, Npath
            Cpath(i,1) = cmplx(normal(),normal(),kind(0d0))*sqrt(0.5d0)*Ampd(i,1)
        end do

        !set H
        do i=1, Nsybl
            do j=1, Npath
                H(i+j-1,i) = Cpath(j,1)
            end do
        end do

        !set HH
        call CAdjoint(H,HH,Nsybl+Npath-1,Nsybl)

        !set HHH
        call CMultiply(HH,H,HHH,Nsybl,Nsybl+Npath-1,Nsybl+Npath-1,Nsybl)

        call CSubstitute(V,HHH,Nsybl,Nsybl)
!        call decomp_zheev(Nsybl,V,Eig)
        call decomp_zheevd(Nsybl,V,Eig)
!        call decomp_zgeev(Nsybl,V,Eig)
!        call decomp_zhpev(Nsybl,V,Eig)

        call sort(Eig,Nsybl)

		do i=1, Nsybl
            output = Eig(1,i)
            OutPutRow = int(output/stride) + 1
            result(i,OutPutRow) = result(i,OutPutRow) + 1.0d0/trial
        end do
	end do
    
    !calculate lambda 1~Nsybl
    do i=1, rank
        lambda(i) = stride*dble(i-1)
    end do
    
    !calculate p(lambda) average lambda
    do i=1, rank
        do j=1, Nsybl
            result(Nsybl+1,i) = result(Nsybl+1,i) + result(j,i)/dble(Nsybl)
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