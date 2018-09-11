program qpsk_ideal
    implicit none

    !declaration
    integer,parameter :: SEbN0=-30
    integer,parameter :: EEbN0=120
    integer KEbN0
    double precision BER
    double precision EbN0

    !initialization
    BER=0.0d0
    EbN0=0.0d0

    !file open
    open (1, file='qpsk_ideal.csv', status='replace')

    !implementation
    do KEbN0=SEbN0, EEbN0
        EbN0 = 10.0d0**(dble(KEbN0/10.0d0)/10.0d0)
        BER = 1.0d0/2.0d0*erfc(sqrt(EbN0))

        write(1,*) dble(kEbN0/10.0d0), ",", BER
    end do

    !file close
    close(1)

end program