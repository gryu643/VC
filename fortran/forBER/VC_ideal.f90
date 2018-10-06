program VC_ideal
    implicit none

    !run time declaration
    integer t1, t2, t_rate, t_max, diff

    !declaration
    integer,parameter :: SEbN0=-10
    integer,parameter :: EEbN0=40
    integer,parameter :: Step=10
    integer i
    integer,parameter :: Nsybl=32
    integer,parameter :: Npath=2
    integer KEbN0
    integer,parameter :: SLambda=0
    integer,parameter :: ELambda=15
    double precision,parameter :: LStep=0.0001
    integer KLambda
    integer SSLambda,EELambda
    double precision BER
    double precision EbN0
    double precision AvLambda(1,8)
    double precision InstantBER
    double precision LambdaEbN0
    double precision ProbLambda
    double precision AvLambdaEbN0
    integer,parameter :: ReadFileRow=150001
    double precision PDF(ReadFileRow,2)
    double precision EbN0dB
    double precision SumPL
    double precision PDF2(ReadFileRow,2)

    !initialization
    BER=0.0d0
    EbN0=0.0d0
    AvLambda=0.0d0
    InstantBER=0.0d0
    LambdaEbN0=0.0d0
    ProbLambda=0.0d0
    AvLambdaEbN0=0.0d0
    PDF=0.0d0
    SumPL=0.0d0

    !average LambdaPDF(each Npath)
    AvLambda(1,1)=0.0d0
    AvLambda(1,2)=0.79d0
    AvLambda(1,4)=0.0d0
    AvLambda(1,8)=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open (1, file='VC_ideal2_0.0001.csv', status='replace')
    open (2, file='evPDFdecomp_s32p2_0.0001.csv', status='old')

    !file read
    do i=1, ReadFileRow
        read(2,*) PDF(i,1), PDF(i,2)
    end do

    SSLambda = nint(dble(SLambda)/LStep)
    EELambda = nint(dble(ELambda)/LStep)

    do KLambda=SSLambda, EELambda
        SumPL = SumPL + LambdaPDF(PDF,ReadFileRow,dble(KLambda)*LStep,Npath)
    end do
!    PDF(:,2) = PDF(:,2)/SumPL

    !implementation
    do KEbN0=SEbN0, EEbN0, Step
        BER=0.0d0
        EbN0 = 10.0d0**(dble(KEbN0)/10.0d0)
        AvLambdaEbN0=0.0d0
!        PDF2(:,1) = PDF(:,1)+EbN0

        do KLambda=SSLambda, EELambda
            LambdaEbN0 = (dble(Klambda)*LStep)*EbN0
            ProbLambda = LambdaPDF(PDF,ReadFileRow,dble(Klambda)*LStep,Npath)/SumPL
            InstantBER = 1.0d0/2.0d0*erfc(sqrt(LambdaEbN0))
            BER = BER + ProbLambda*InstantBER
            AvLambdaEbN0 = AvLambdaEbN0 + LambdaEbN0*ProbLambda
        end do

        EbN0dB = 10.0d0*dlog10(AvLambdaEbN0)
        print *, "AvEbN0=", EbN0dB
        print *, "BER=", BER
        write(1,*) EbN0dB, ",", BER
    end do

    !file close
    close(1)
    close(2)

    !time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)

contains
    function LambdaPDF(PDF,Row,x,path)
        !-- declaration
        double precision x
        integer path
        integer i
        double precision LambdaPDF
        double precision y
        double precision PDF(Row,2)
        integer Row

        !-- initialization
        LambdaPDF=0.0d0
        y=0.0d0
            
        do i=1, Row
            if((PDF(i,1)<=x).and.(x<PDF(i+1,1))) then
                y = PDF(i,2)
                exit
            endif
        end do

        LambdaPDF = y
    end function LambdaPDF

end program VC_ideal
