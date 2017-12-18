
! Free form main program

program PPL_Soturon
PARAMETER(lcar=1024,lbit=512,lsym=128,ldat=8192,lndr=512,lpath=12,lrow=30)
COMPLEX C(lrow,lrow),BRAMDA(lrow,lrow),C1(lrow,lrow),C2(lrow,lrow)
COMPLEX C3(lrow,lrow),C4(lrow,lrow),C5(lrow,lrow)
COMPLEX D(lrow,1),D1(lrow,1),D3(lrow,1),D4(lrow,1),D5(lrow,1)
COMPLEX D6(lrow,1),D7(lrow,1),D8(lrow,1),D9(lrow,1),E(1,lrow),G(lrow-1,1)
COMPLEX CLAMDA(lrow,lrow),P(lrow,lrow),S(lrow,lrow),Q(500,1)
COMPLEX F(lrow+lpath,1),HC(lrow+lpath,lrow),HC2(lrow+lpath,lrow)
COMPLEX H(lrow+lpath,lrow),HE(lrow+2*lpath,lrow+lpath)
COMPLEX HEC(lrow+2*lpath,lrow),HH(lrow,lrow+lpath),HHH(lrow,lrow)
COMPLEX N(1,1),CH(lrow,lrow),CTMP(lrow,lrow)
INTEGER ROWB,U
REAL    D2(lrow),M1,M2,YYY
character TMP*(39),FNAME*(30)

read(*,*)TMP,ROWB !ƒVƒ“ƒ{ƒ‹”
read(*,*)TMP,U !ƒpƒX”
read(*,*)TMP,FNAME !o—Íƒtƒ@ƒCƒ‹–¼
OPEN(24,FILE=FNAME,STATUS='UNKNOWN',ACCESS='SEQUENTIAL')

H(:,:)=cmplx(0.,0.)
HE(:,:)=cmplx(0.,0.)
HHH(:,:)=cmplx(0.,0.)
HC(:,:)=cmplx(0.,0.)
HC2(:,:)=cmplx(0.,0.)
HEC(:,:)=cmplx(0.,0.)
BRAMDA(:,:)=cmplx(0.,0.)
CLAMDA(:,:)=cmplx(0.,0.)
P(:,:)=cmplx(0.,0.)
S(:,:)=cmplx(0.,0.)
Q(:,:)=cmplx(0.,0.)
CH(:,:)=cmplx(0.,0.)
N(:,:)=cmplx(0.,0.)
C(:,:)=cmplx(0.,0.)
C1(:,:)=cmplx(0.,0.)
C2(:,:)=cmplx(0.,0.)
C3(:,:)=cmplx(0.,0.)
C4(:,:)=cmplx(0.,0.)
C5(:,:)=cmplx(0.,0.)
D(:,:)=cmplx(0.,0.)
D1(:,:)=cmplx(0.,0.)
D3(:,:)=cmplx(0.,0.)
D4(:,:)=cmplx(0.,0.)
D5(:,:)=cmplx(0.,0.)
D6(:,:)=cmplx(0.,0.)
D7(:,:)=cmplx(0.,0.)
D8(:,:)=cmplx(0.,0.)
D9(:,:)=cmplx(0.,0.)
E(:,:)=cmplx(0.,0.)
G(:,:)=cmplx(0.,0.)
D2(:)=0.




!“`”À˜Hs—ñH‚ðÝ’è(,ROWB)
do j=0,U-1
	do i=1,ROWB
		H(i+j,i)=cmplx(0.1+0.1*j,0.2+0.1*j)
	end do
end do


!“`”À˜Hs—ñH‚ðŠg’£‚µ‚½HE‚ðl‚¦‚é
do j=0,U-1
	DO i=1,ROWB+U-1
		HE(i+j,i)=cmplx(0.1+0.1*j,0.2+0.1*j)
	end do
end do


!”º“`”À˜Hs—ñHH‚ðÝ’è
CALL CADJOINT(H,HH,lrow+lpath,lrow)
!HH~H
CALL CAMULTIPLY(HH,H,HHH,lrow,lrow+lpath,lrow+lpath,lrow)


!”CˆÓ“`‘—ƒxƒNƒgƒ‹‚ðÝ’è
c(:,:)=cmplx(0.,0.)

do j=1,ROWB
	do i=1,ROWB
		c(i,j)=cmplx(1.,0.)
	end do
end do


do l=1,200

c(:,:)=cmplx(0.,0.)

do j=1,ROWB
	do i=1,ROWB
		c(i,j)=cmplx(1.,0.)
	end do
end do


	M1=0
	M2=0
	do k=1,l
		
		!ŽŸƒ‹[ƒv‚Å‚ÌŒÅ—L’lŽZo‚Ì‚½‚ßAC1‚É‘Oƒ‹[ƒv‚ÌŒÅ—LƒxƒNƒgƒ‹ŒQ‚ðŠi”[
		C1(:,:)=C(:,:)
		!ƒƒ“ƒ‹[ƒvŠÖ”‚ð—ñƒxƒNƒgƒ‹ŒQ‚É“K‰ž‚³‚¹‚é
		CALL CAMULTIPLY(H,C,HC,lrow+lpath,lrow,lrow,lrow)
	
		
		!HC‚Éˆ—I‚ð‰Á‚¦‚é
		do j=1,ROWB
			do i=1,ROWB+U-1
				HC2(i,j)=conjg(HC(ROWB+U-i,j))
			end do
		end do
		
		!HE’Ê‰ß
		CALL CAMULTIPLY(HE,HC2,HEC,lrow+2*lpath,lrow+lpath,lrow+lpath,lrow)
				
		!HEC‚Éˆ—J‚ð‰Á‚¦‚é
		do j=1,ROWB
			do i=1,ROWB
				C(i,j)=conjg(HEC(ROWB+U-i,j))
			end do
		end do

!		do i=1, ROWB
!			print *, l, C(i,2)
!		end do

		!—ñƒxƒNƒgƒ‹ŒQ‚ÌŒÅ—L’l‚ð‚»‚ê‚¼‚êŽZo
		do j=1,ROWB
	
			do i=1,ROWB
				D(i,1)=C(i,j) !–³ü‹@ŠÖ‰•œŒã‚Ì—ñƒxƒNƒgƒ‹‚ðŠi”[
				D1(i,1)=C1(i,j) !–³ü‹@ŠÖ‰•œ‘O‚Ì—ñƒxƒNƒgƒ‹‚ðŠi”[
				D2(i)=abs(real(D(i,1))/real(D1(i,1))) !‰•œŒã‚ÌƒxƒNƒgƒ‹‚ð‘O‚Ì‚ÅŠ„‚Á‚ÄAŒÅ—L’l‚Æ‚·‚é
			end do
			
			blamda=D2(1)
			do i=2,ROWB
				if(blamda.GT.D2(i)) blamda=D2(i) !ã‚Å‹‚ß‚½Šes‚ÌŒÅ—L’l‚Ì’†‚ÅÅ¬‚Ì‚à‚Ì‚ðblamda‚ÉŠi”[
			end do
			bramda(:,j)=cmplx(blamda,0.) !ŒÅ—L’l‚ð•¡‘f”Œ`‚Åbramda‚Ì—ñ‚ÉˆêŠ‡Ši”[
		end do
	
!		do i=1, ROWB
!			print *, l, bramda(1,i)
!		end do	
				
		!·•ª“±o
		C2(:,:)=cmplx(0.,0.)
		do j=1,ROWB-1 !j=1~4
			do i=1,ROWB !i=1~5
				D3(i,1)=C1(i,j) !–³ü‹@ŠÖ‰•œ‘O‚ÌŒÅ—LƒxƒNƒgƒ‹
				D4(i,1)=C1(i,j+1) !D3‚ÌŽŸ‚Ì—ñ‚ÌŒÅ—LƒxƒNƒgƒ‹
			end do
			
			CALL CADJOINT(D3,E,lrow,1)
			
			do i=1,ROWB
				D5(i,1)=cmplx(real(bramda(i,j))*real(D3(i,1)),real(bramda(i,j))*aimag(D3(i,1)))
			end do
			
			CALL CAMULTIPLY(D5,E,C3,lrow,1,1,lrow)
			
			CALL CASUM(C2,C3,C2,lrow,lrow)
			
			CALL CAMULTIPLY(C2,D4,D6,lrow,lrow,lrow,1)
			
			do i=1,ROWB
				C4(i,j)=D6(i,1)
			end do
		end do
		
		do j=1,ROWB-1
			do i=1,ROWB
				C5(i,j+1)=C4(i,j)
			end do
		end do
		
		CALL CASUB(C,C5,C,lrow,lrow)

!		do i=1, ROWB
!			print *, l, C5(i,4)
!		end do
		
		!³‹K‰»(‚±‚±‚Å‚Í—ñƒxƒNƒgƒ‹ŒQ“à‚Ì—ñƒxƒNƒgƒ‹‚ðˆê‚Â‚¸‚Â³‹K‰»‚µ‚Ä‚¢‚é)
		do j=1,ROWB
			do i=1,ROWB
				D7(i,1)=C(i,j)
			end do

			CALL CNOM(D7,lrow,1)
		
			do i=1,ROWB
				C(i,j)=D7(i,1)
			end do
		end do

!		do i=1, ROWB
!			print *, l, C(i,4)
!		end do		

	END DO

!	do i=1, ROWB
!		print *, l, C(i,1)
!	end do

	!ŒÅ—LƒxƒNƒgƒ‹‚©Šm”F(“àÏ=0)
	G(1,1)=cmplx(0.0,0.0)
	do i=1,ROWB
		do j=i+1, ROWB
			do k=1,ROWB
				D8(k,1)=C(k,i)
				D9(k,1)=C(k,j)
			end do
			CALL CHN(D8,D9,N,lrow,1)
			G(1,1)=N(1,1)
		end do
	end do

	call CNORM(G,YYY,1,1)
	print *,l,",",YYY
	
	
	!ŒŸØ
	!ƒXƒyƒNƒgƒ‹’è—
	do i=1,ROWB
		CLAMDA(i,i)=bramda(1,i)
	end do
	CALL CAMULTIPLY(C,CLAMDA,CTMP,lrow,lrow,lrow,lrow)
	CALL CADJOINT(C,CH,lrow,lrow)
	CALL CAMULTIPLY(CTMP,CH,P,lrow,lrow,lrow,lrow)
	
	do j=1,ROWB
		do i=1,ROWB
			P(i,j)=cmplx(abs(real(P(i,j))),aimag(P(i,j)))
		end do
	end do
	
	!s—ñ‚Ì·‚Ì“ñæ˜a
	CALL CASUB(HHH,P,S,lrow,lrow)

	do j=1,ROWB
		do i=1,ROWB
			M1=M1+real(S(i,j))**2+aimag(S(i,j))**2
			M2=M2+real(HHH(i,j))**2+aimag(HHH(i,j))**2
		end do
	end do
	Q(l,1)=M1/M2
	
END DO

do i=1,200
	write(24,*)i,real(Q(i,1))
end do



end program

SUBROUTINE CADJOINT(A,AH,ldat,lcar)
COMPLEX A(ldat,lcar),AH(lcar,ldat)

do i=1,ldat
	do j=1,lcar
	AH(j,i)=conjg(A(i,j))
	end do
end do
END

SUBROUTINE CAMULTIPLY(A,B,C,ldat,lcar,ldat2,lcar2)
COMPLEX A(ldat,lcar),B(ldat2,lcar2),C(ldat,lcar2)
IF(lcar.NE.ldat2)then
	write(*,*)'error in array multiplication'
	STOP
end if
 C(:,:)=cmplx(0.,0.)
do i=1,ldat
	do j=1,lcar2
		do k=1,ldat2
			C(i,j)=C(i,j)+A(i,k)*B(k,j)
		end do
	end do
end do

END

SUBROUTINE CASUM(A,B,C,ldat,lcar)
complex A(ldat,lcar),B(ldat,lcar),C(ldat,lcar)

do i=1,ldat
	do j=1,lcar
		C(i,j)=A(i,j)+B(i,j)
	end do
end do

END 

SUBROUTINE CASUB(A,B,C,ldat,lcar)
complex A(ldat,lcar),B(ldat,lcar),C(ldat,lcar)

do i=1,ldat
	do j=1,lcar
		C(i,j)=A(i,j)-B(i,j)
	end do
end do

END
	
		
SUBROUTINE CNORM(A,XOUT,ldat,lcar)
complex A(ldat,lcar)
real XOUT

XOUT=0.
do i=1,ldat
	XOUT=XOUT+real(A(i,1))**2+aimag(A(i,1))**2
end do

XOUT=sqrt(XOUT)

END

SUBROUTINE CNOM(A,ldat,lcar)
complex A(ldat,lcar)
real XIN

CALL CNORM(A,XIN,ldat,lcar)


do i=1,ldat
	A(i,1)=A(i,1)/XIN
end do

END

SUBROUTINE CHN(A,B,C,ldat,lcar)
COMPLEX A(ldat,lcar),B(ldat,lcar),C(lcar,lcar),X(lcar,ldat)

do i=1,ldat
	X(1,i)=conjg(A(i,1))
end do

CALL CAMULTIPLY(X,B,C,lcar,ldat,ldat,lcar)

END