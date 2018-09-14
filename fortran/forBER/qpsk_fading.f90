program qpsk_fading
    use CALmod
    implicit none

    !run time declaration
    integer t1, t2, t_rate, t_max, diff

    !declaration
    integer,parameter :: SEbN0=-0
    integer,parameter :: EEbN0=40
    integer,parameter :: Step=10
    integer,parameter :: Nloop=100000

    integer i,j
    double precision Ps
    double precision Pn
    integer Collect
    integer False
    integer loop
    integer KEbN0
    complex(kind(0d0)) S
    integer TdatI
    integer TdatQ
    double precision Pow
    complex(kind(0d0)) Noise
    complex(kind(0d0)) Rayl
    integer RdatI
    integer RdatQ
    double precision EbN0
    double precision BER

    !initialize
    Ps=0.0d0
    Pn=0.0d0
    Collect=0
    False=0
    loop=0
    KEbN0=0
    S=0.0d0
    TdatI=0
    TdatQ=0
    Pow=0.0d0
    Noise=(0.0d0,0.0d0)
    Rayl=(0.0d0,0.0d0)
    RdatI=0
    RdatQ=0
    EbN0=0.0d0
    BER=0.0d0

    !time measurement start
    call system_clock(t1)

    !file open
    open (1, file='qpsk_fading.csv', status='replace')

    !implimentation part
    do KEbN0=SEbN0,EEbN0, Step !Ps/Pn loop
        Ps = 0.0d0
        Pn = 0.0d0
        Collect = 0
        False = 0

        do loop=1, Nloop !Monte calro loop
            !set information symbol
            S = cmplx(nint(rand())*2.0d0-1.0d0,nint(rand())*2.0d0-1.0d0,kind(0d0)) !-1or1
            TdatI = (real(S)+1)/2 !0or1
            TdatQ = (aimag(S)+1)/2

            !signal power = 1
            Pow = 0.0d0
            Pow = abs(S)**2.0d0
            S = S / sqrt(Pow)

            !generate Noise
            Noise = cmplx(normal(),normal(),kind(0d0))

            !generate rayleigh
            Rayl = cmplx(normal(),normal(),kind(0d0))/sqrt(2.0d0)/sqrt(8.0d0)

            !EbN0
            Noise = Noise * sqrt(1.0d0/(10.0d0**(KEbN0/10.0d0))/2.0d0)/sqrt(2.0d0)

            !multiply rayleigh
            S = Rayl * S

            !calculate Ps and Pn
            !Noiseが統計的に電力１になるので、試行回数分電力を足し合わせる
            Ps = Ps + abs(S)**2.0d0
            Pn = Pn + abs(Noise)**2.0d0

            !add noise
            S = S + Noise

            RdatI = 0
            RdatQ = 0
            !Demodulation
            if(real(S)>0.0) then
                RdatI = 1
            elseif(real(S)<0.0) then
                RdatI = 0
            endif
            if(aimag(S)>0.0) then
                RdatQ = 1
            elseif(aimag(S)<0.0) then
                RdatQ = 0
            endif

            !Bit Error Rate
            if(RdatI==TdatI) then
                Collect = Collect + 1
            elseif(RdatI/=TdatI) then
                False = False + 1
            endif
            if(RdatQ==TdatQ) then
                Collect = Collect + 1
            elseif(RdatQ/=TdatQ) then
                False = False + 1
            endif
        end do
        
        EbN0 = 10.0d0*dlog10(Ps/Pn/2.0d0) !1/2*SNR=EbN0 (QPSK)
        BER = dble(False) / (dble(Collect) + dble(False))
        if(BER>0.0) then
            write(1,*) EbN0, ',', BER
            print *, 'EbN0=', EbN0
            print *, 'BER=', BER
        endif
    end do

    !file close
    close(1)

    !time measurement end
    call system_clock(t2, t_rate, t_max)
    print *, 'Elapsed time is...', (t2-t1)/dble(t_rate)
end program