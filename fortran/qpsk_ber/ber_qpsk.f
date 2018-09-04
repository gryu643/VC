      PROGRAM QAM BER CAL               

      REAL PES(300),A

      DO I=-30,120
        A=10.**(float(I)/100.)

        A=SQRT(A)                        
        PES(I)=1./2.*ERFC(A) 

      WRITE(10,*)I/10.,'	',PES(I)
      end do

      STOP
      END


      REAL FUNCTION ERFC(X)
      PI=4.*ATAN(1.)
      X2=ABS(X)*SQRT(2.)
      IF(X2 .LE. 2.5) THEN
         X3=1.
         X4=3.
         X5=2.
         X6=1.
         X7=1.
         X8=0.
         X9=1.
10       X10=-x2**(2.*X3)/(X4*X5*X6*X7)
         X8=X8+X10
         X11=ABS((1.+X8-X9)/(1. +X8))
         X9=1.+X8
         IF(X11 .LE. 1.E-06) THEN   !
           IF(X .GE. 0.) THEN
              ERFC=1.-2.*X2*(1.+X8)/SQRT(2.*PI)
           ELSE 
              ERFC=1.+2.*X2*(1.+X8)/SQRT(2.*PI)
           END IF
           RETURN
         END IF
         X4=X4+2.
         X5=-2.*X5
         X6=X6*X7
         X7=X7+1.
         X3=X3+1.
         GO TO 10
       END IF
       X3=10.
       X4=1.
       X5=1.
20     X3=X3+1.
       X6=X2+X3
30     X6=X2+(X3-X4)/X6
       IF((X3-X4) .NE. 1.) THEN
           X4=X4+1.
           GO TO 30
       END IF
       X6=1./X6
       X7=ABS((X6-X5)/X6)
       X5=X6
       IF(X7 .LE. 1.E-6) THEN
          IF(X. GE. 0.) THEN
            ERFC=2.*EXP(-X2*X2/2.)*X6/SQRT(2.*PI)
          ELSE
            ERFC=2.-2.*EXP(-X2*X2/2.)*X6/SQRT(2.*PI)
          END IF
          RETURN
        END IF
        X4=1.
        GO TO 20
        END


