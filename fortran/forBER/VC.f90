program VC
    use PPLmod
    use CALmod
    implicit none

    !declaration
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=8
    integer,parameter :: SEbN0=-10
    integer,parameter :: EEbN0=40
    integer,parameter :: Step=5
    integer,parameter :: Nloop=10000
    integer,parameter :: PPLloop=500

    integer i
    double precision Ampd(Npath,1)
    double precision Psign
    double precision Pwgn
    integer Collect
    integer False
    integer loop
    integer EbN0
    complex(kind(0d0)) Cpath(Npath,1)
    complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
    complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
    complex(kind(0d0)) X(Nsybl,Nsybl)

    !initialize
    Ampd(:,:)=0.0
    Psign=0.0
    Pwgn=0.0
    Collect=0
    False=0
    loop=0
    EbN0=0
    Cpath(:,:)=(0.0,0.0)
    H(:,:)=(0.0,0.0)
    HE(:,:)=(0.0,0.0)
    X(:,:)=(0.0,0.0)


    !implimentation part

    !channel gain parameter
    do i=1, Npath
        Ampd(i,1) = sqrt(1/Npath) !Equal Gain
        !Ampd(i) = sqrt(1/(2^(i-1))) !Exp. atten.
    end do

    do EbN0=SEbN0, EEbN0, Step
        Psign = 0.0
        Pwgn = 0.0
        Collect = 0
        False = 0

        do loop=1, Nloop
            do i=1, Npath
                Cpath(i) = cmplx(normal(),normal(),kind(0d0))/sqrt(2)*Ampd(i)
            end do
        end do 

        do i=1, Nsybl
            do j=1, Npath
                H(i+j-1,i) = Cpath(j)
            end do
        end do

        do i=1, Nsybl+Npath-1
            do j=1, Npath
                HE(i+j-1,i) = Cpath(j)
            end do
        end do

        do i=1, Nsybl
            do j=1, Nsybl
                X(i,j) = cmplx(1.0, 0.0, kind(0d0))
            end do
        end do
    
    end do

    !file open
    open (1, file='VC(Nsybl=32,Npath=8).csv', status='replace')


    !file close
    close(18)
end program