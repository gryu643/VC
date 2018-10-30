program eigenvaluePDF
    use CALmod
    implicit none

    !run time declaration
    integer t1, t2, t_rate, t_max, diff

    !declatation
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=8
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

    double precision,parameter :: stride=0.0001d0
    double precision output
    double precision,allocatable :: result(:,:)
    double precision,allocatable :: lambda(:)
    integer i,j,loop
    double precision,parameter :: start=0.0d0
    double precision,parameter :: end=25.0d0
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
    rank=nint((end-start)/stride)

    !allocate
    allocate(result(Nsybl+1,rank))
    allocate(lambda(rank))

    !allocate initialize
    result=0.0d0
    lambda=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open(7,file='ev(s32p8)L1e-4.csv', status='replace')
    open(8,file='ev(s32p8)L2e-4.csv', status='replace')
    open(9,file='ev(s32p8)L3e-4.csv', status='replace')
    open(10,file='ev(s32p8)L4e-4.csv', status='replace')
    open(11,file='ev(s32p8)L5e-4.csv', status='replace')
    open(12,file='ev(s32p8)L6e-4.csv', status='replace')
    open(13,file='ev(s32p8)L7e-4.csv', status='replace')
    open(14,file='ev(s32p8)L8e-4.csv', status='replace')
    open(15,file='ev(s32p8)L9e-4.csv', status='replace')
    open(16,file='ev(s32p8)L10e-4.csv', status='replace')
    open(17,file='ev(s32p8)L11e-4.csv', status='replace')
    open(18,file='ev(s32p8)L12e-4.csv', status='replace')
    open(19,file='ev(s32p8)L13e-4.csv', status='replace')
    open(20,file='ev(s32p8)L14e-4.csv', status='replace')
    open(21,file='ev(s32p8)L15e-4.csv', status='replace')
    open(22,file='ev(s32p8)L16e-4.csv', status='replace')
    open(23,file='ev(s32p8)L17e-4.csv', status='replace')
    open(24,file='ev(s32p8)L18e-4.csv', status='replace')
    open(25,file='ev(s32p8)L19e-4.csv', status='replace')
    open(26,file='ev(s32p8)L20e-4.csv', status='replace')
    open(27,file='ev(s32p8)L21e-4.csv', status='replace')
    open(28,file='ev(s32p8)L22e-4.csv', status='replace')
    open(29,file='ev(s32p8)L23e-4.csv', status='replace')
    open(30,file='ev(s32p8)L24e-4.csv', status='replace')
    open(31,file='ev(s32p8)L25e-4.csv', status='replace')
    open(32,file='ev(s32p8)L26e-4.csv', status='replace')
    open(33,file='ev(s32p8)L27e-4.csv', status='replace')
    open(34,file='ev(s32p8)L28e-4.csv', status='replace')
    open(35,file='ev(s32p8)L29e-4.csv', status='replace')
    open(36,file='ev(s32p8)L30e-4.csv', status='replace')
    open(37,file='ev(s32p8)L31e-4.csv', status='replace')
    open(38,file='ev(s32p8)L32e-4.csv', status='replace')
    open(39,file='ev(s32p8)e-4.csv', status='replace')

    !implementation
    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1.0d0/dble(Npath)) !Equal Gain
        !Ampd(i) = sqrt(1/(2**(i-1))) !Exp. atten.
    end do


	do loop=1, trial
        if(mod(loop,int(trial/10.0d0))==0) then
            print *, loop
        endif
        
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

        !decomp
        if(Npath==1) then
            do i=1, Nsybl
                Eig(1,i) = HHH(i,i)
                V(i,i) = cmplx(1.0d0,0.0d0,kind(0d0))
            end do
        else
            call CSubstitute(V,HHH,Nsybl,Nsybl)
    !        call decomp_zheev(Nsybl,V,Eig)
            call decomp_zheevd(Nsybl,V,Eig)
    !        call decomp_zgeev(Nsybl,V,Eig)
    !        call decomp_zhpev(Nsybl,V,Eig)
        endif

        call sort(Eig,Nsybl)

		do i=1, Nsybl
            output = Eig(1,i)
            OutPutRow = int(output/stride) + 1
            result(i,OutPutRow) = result(i,OutPutRow) + 1.0d0/trial
        end do
	end do
    
    !calculate lambda 1~Nsybl
    do i=1, rank
        lambda(i) = stride*dble(i)
    end do
    
    !calculate p(lambda) average lambda
    do i=1, rank
        do j=1, Nsybl
            result(Nsybl+1,i) = result(Nsybl+1,i) + result(j,i)/dble(Nsybl)
        end do
    end do

    do i=1, Nsybl+1
        do j=1, rank
            write(i+6,*) lambda(j), ',', result(i,j)
        end do
    end do

    !file close
    do i=7,Nsybl+7
        close(i)
    end do
    
    !time measurement end
    call system_clock(t2,t_rate,t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program