program VC_ideal
    implicit none

    !run time declaration
    integer t1, t2, t_rate, t_max, diff

    !declaration
    integer,parameter :: SEbN0=-10
    integer,parameter :: EEbN0=40
    integer,parameter :: Step=10
    integer i,j
    integer,parameter :: Nsybl=32
    integer KEbN0
    integer,parameter :: SLambda=0
    integer,parameter :: ELambda=20
    !入力する確率密度分布の、横軸刻み幅に合わせる
    double precision,parameter :: LStep=0.0001
    double precision,parameter :: LLStep=0.0001

    integer KLambda
    integer SSLambda,EELambda
    double precision BER
    double precision EbN0
    double precision InstantBER
    double precision LambdaEbN0
    double precision ProbLambda
    double precision AvLambdaEbN0
    !入力する確率密度分布の、lambda=15までの行数
    integer,parameter :: ReadFileRow=int(ELambda/Lstep)
    double precision PDF(ReadFileRow,2)
    double precision PDFother(ReadFileRow,Nsybl+1)
    double precision EbN0dB
    double precision SumPL
    double precision SumPLother(1,Nsybl)
    double precision AvLambda(1,Nsybl)
    double precision AL

    !initialization
    BER=0.0d0
    EbN0=0.0d0
    InstantBER=0.0d0
    LambdaEbN0=0.0d0
    ProbLambda=0.0d0
    AvLambdaEbN0=0.0d0
    PDF=0.0d0
    PDFother=0.0d0
    SumPL=0.0d0
    SumPLother=0.0d0
    AvLambda=0.0d0
    AL=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open(7,file='ev(s32p8)1Le-4.csv', status='old')
    open(8,file='ev(s32p8)2Le-4.csv', status='old')
    open(9,file='ev(s32p8)3Le-4.csv', status='old')
    open(10,file='ev(s32p8)4Le-4.csv', status='old')
    open(11,file='ev(s32p8)5Le-4.csv', status='old')
    open(12,file='ev(s32p8)6Le-4.csv', status='old')
    open(13,file='ev(s32p8)7Le-4.csv', status='old')
    open(14,file='ev(s32p8)8Le-4.csv', status='old')
    open(15,file='ev(s32p8)9Le-4.csv', status='old')
    open(16,file='ev(s32p8)10Le-4.csv', status='old')
    open(17,file='ev(s32p8)11Le-4.csv', status='old')
    open(18,file='ev(s32p8)12Le-4.csv', status='old')
    open(19,file='ev(s32p8)13Le-4.csv', status='old')
    open(20,file='ev(s32p8)14Le-4.csv', status='old')
    open(21,file='ev(s32p8)15Le-4.csv', status='old')
    open(22,file='ev(s32p8)16Le-4.csv', status='old')
    open(23,file='ev(s32p8)17Le-4.csv', status='old')
    open(24,file='ev(s32p8)18Le-4.csv', status='old')
    open(25,file='ev(s32p8)19Le-4.csv', status='old')
    open(26,file='ev(s32p8)20Le-4.csv', status='old')
    open(27,file='ev(s32p8)21Le-4.csv', status='old')
    open(28,file='ev(s32p8)22Le-4.csv', status='old')
    open(29,file='ev(s32p8)23Le-4.csv', status='old')
    open(30,file='ev(s32p8)24Le-4.csv', status='old')
    open(31,file='ev(s32p8)25Le-4.csv', status='old')
    open(32,file='ev(s32p8)26Le-4.csv', status='old')
    open(33,file='ev(s32p8)27Le-4.csv', status='old')
    open(34,file='ev(s32p8)28Le-4.csv', status='old')
    open(35,file='ev(s32p8)29Le-4.csv', status='old')
    open(36,file='ev(s32p8)30Le-4.csv', status='old')
    open(37,file='ev(s32p8)31Le-4.csv', status='old')
    open(38,file='ev(s32p8)32Le-4.csv', status='old')
    open(39,file='ev(s32p8)e-4.csv', status='old')
    open (1, file='VC_ideal(s32p8)Be-4.csv', status='replace')

    !file read
    do i=1, ReadFileRow
        read(39,*) PDF(i,1), PDF(i,2)
    end do
    do i=1, Nsybl
        do j=1, ReadFileRow
            read(i+6,*) PDFother(j,1), PDFother(j,i+1)
        end do
    end do
    !-- implementation
    SSLambda = nint(dble(SLambda)/LLStep)+1
    EELambda = nint(dble(ELambda)/LLStep)

    do KLambda=SSLambda, EELambda
        SumPL = SumPL + LambdaPDF(PDF,ReadFileRow,dble(KLambda)*LLStep,LStep)
    end do

    !calculate average lambda
    do i=1, Nsybl
        do j=1, ReadFileRow
            AvLambda(1,i) = AvLambda(1,i) + PDFother(j,1)*PDFother(j,i+1)
        end do
    end do
    do i=1, Nsybl
        AL = AL + AvLambda(1,i)/Nsybl
    end do

    do KEbN0=SEbN0, EEbN0, Step
        BER=0.0d0
        EbN0 = 10.0d0**(dble(KEbN0)/10.0d0)
        AvLambdaEbN0=0.0d0

        do KLambda=SSLambda, EELambda
            LambdaEbN0 = (dble(Klambda)*LLStep)*EbN0
            ProbLambda = LambdaPDF(PDF,ReadFileRow,dble(Klambda)*LLStep,LStep)/SumPL
            InstantBER = 1.0d0/2.0d0*erfc(sqrt(LambdaEbN0))
            BER = BER + ProbLambda*InstantBER
!            AvLambdaEbN0 = AvLambdaEbN0 + LambdaEbN0*ProbLambda
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
    function LambdaPDF(PDF,Row,x,LStep)
        !-- declaration
        double precision x
        double precision LambdaPDF
        double precision y2,y1
        integer x2,x1
        double precision PDF(Row,2)
        integer Row
        double precision LStep
        double precision Slope

        !-- initialization
        LambdaPDF=0.0d0
            
        !-- implementation
        x1 = int(x/LStep)+1
        x2 = x1+1
        y1 = PDF(x1,2)
        y2 = PDF(x2,2)
        Slope = (y2-y1)/LStep

        LambdaPDF = y1 + Slope*(x-(PDF(x1,1)-LStep*0.5d0))
    end function LambdaPDF

end program VC_ideal
