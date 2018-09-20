program qpsk_ideal_MRdiv
    implicit none

    !declaration
    integer,parameter :: SEbN0=-100
    integer,parameter :: EEbN0=400
    integer KEbN0
    double precision BER
    double precision EbN0

    !initialization
    BER=0.0d0
    EbN0=0.0d0

    !file open
    open (1, file='qpsk_ideal_MRdiv.csv', status='replace')

    !implementation
    do KEbN0=SEbN0, EEbN0
        EbN0 = 10.0d0**(dble(KEbN0/10.0d0)/10.0d0)
        BER = 1.0d0/2.0d0 - 1.0d0/2.0d0*(1.0d0+3.0d0/EbN0)/(1.0d0+2.0d0/EbN0)**1.5d0

        write(1,*) dble(kEbN0/10.0d0+1.6825), ",", BER
    end do

    !file close
    close(1)

end program