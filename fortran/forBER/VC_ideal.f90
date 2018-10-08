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
    integer KEbN0
    integer,parameter :: SLambda=0
    integer,parameter :: ELambda=20
    !入力する確率密度分布の、横軸刻み幅に合わせる
    double precision,parameter :: LStep=0.000001

    integer KLambda
    integer SSLambda,EELambda
    double precision BER
    double precision EbN0
    double precision InstantBER
    double precision LambdaEbN0
    double precision ProbLambda
    double precision AvLambdaEbN0
    !入力する確率密度分布の、lambda=15までの行数
    integer,parameter :: ReadFileRow=int(ELambda/Lstep)+1
    double precision PDF(ReadFileRow,2)
    double precision EbN0dB
    double precision SumPL

    !initialization
    BER=0.0d0
    EbN0=0.0d0
    InstantBER=0.0d0
    LambdaEbN0=0.0d0
    ProbLambda=0.0d0
    AvLambdaEbN0=0.0d0
    PDF=0.0d0
    SumPL=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open (1, file='VC_ideal(s32p2)e-6.csv', status='replace')
    open (2, file='ev(s32p2)Ae-6.csv', status='old')

    !file read
    do i=1, ReadFileRow
        read(2,*) PDF(i,1), PDF(i,2)
    end do

    SSLambda = nint(dble(SLambda)/LStep)
    EELambda = nint(dble(ELambda)/LStep)

    do KLambda=SSLambda, EELambda
        SumPL = SumPL + LambdaPDF(PDF,ReadFileRow,dble(KLambda)*LStep,LStep)
    end do

    !implementation
    do KEbN0=SEbN0, EEbN0, Step
        BER=0.0d0
        EbN0 = 10.0d0**(dble(KEbN0)/10.0d0)
        AvLambdaEbN0=0.0d0

        do KLambda=SSLambda, EELambda
            LambdaEbN0 = (dble(Klambda)*LStep)*EbN0
            ProbLambda = LambdaPDF(PDF,ReadFileRow,dble(Klambda)*LStep,LStep)/SumPL
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
    function LambdaPDF(PDF,Row,x,Step)
        !-- declaration
        double precision x
        integer i
        double precision LambdaPDF
        double precision y
        double precision PDF(Row,2)
        integer Row
        double precision Step

        !-- initialization
        LambdaPDF=0.0d0
        y=0.0d0
            
        i = int(x/Step) + 1
        y = PDF(i,2)

        LambdaPDF = y
    end function LambdaPDF

end program VC_ideal
