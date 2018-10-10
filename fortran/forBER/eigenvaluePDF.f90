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

    double precision,parameter :: stride=0.00001d0
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
    open(39,file='ev(s32p8)Be-5.csv', status='replace')

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
        lambda(i) = stride*dble(i)
    end do
    
    !calculate p(lambda) average lambda
    do i=1, rank
        do j=1, Nsybl
            result(Nsybl+1,i) = result(Nsybl+1,i) + result(j,i)/dble(Nsybl)
        end do
    end do

    do j=1, rank
        write(39,*) lambda(j), ',', result(33,j)
    end do

    !file close
    close(39)
    
    !time measurement end
    call system_clock(t2,t_rate,t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program