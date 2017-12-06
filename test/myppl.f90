program ppl_bulk
	implicit none

	integer,parameter :: H_ROW=1024,H_COL=1024,X_ROW=1024,X_COL=1024,H_PATH=64
	integer i,j,l,m
	integer SYMBL,PATH
	character(10) OUT_FNAME
	character(10) TMP
	complex :: X(X_ROW,X_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex :: H(H_ROW,H_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex :: HE(H_ROW+H_PATH,H_COL)=(0.0,0.0) !(SYMBL+PATH-1,SYMBL)
	complex :: HH(H_ROW,H_COL+H_PATH)=(0.0,0.0) !(SYMBL,SYMBL+PATH-1)
	complex :: HHH(H_ROW,H_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex :: HEX(H_ROW+H_PATH,X_COL)=(0.0,0.0) !(SYMBL+PATH-1,SYMBL)

	!シンボル数、パ数数、出力ファイル名を読み込む
	read(5,*) TMP,SYMBL
	read(5,*) TMP,PATH
	read(5,*) TMP,OUT_FNAME

	!伝搬路行列Hの設定
	do i=1, SYMBL
		do j=1, PATH
			if((i+j-1).gt.SYMBL) exit
			H(i+j-1,i) = cmplx(0.1*j, 0.1+0.1*j)
		end do
	end do

	!伝搬路行列Hを拡張したHEを設定
	do i=1, SYMBL
		do j=1, PATH
			HE(i+j-1,i) = cmplx(0.1*j, 0.1+0.1*j)
		end do
	end do

	!行列Hの随伴行列HHの設定
	call CmplxAdjoint(H,HH,SYMBL,SYMBL)

	!合成チャネル行列HHHの設定
	call CmplxMultiply(HH,H,HHH,SYMBL,SYMBL,SYMBL,SYMBL)

	do l=1, 200
		!任意伝送ベクトルの設定
		do i=1, SYMBL
			do j=1, SYMBL
				X(i,j) = cmplx(1.0, 0.0)
			end do
		end do

		do m=1, l
			!伝搬路HE通過
			call CmplxMultiply(HE,X,HEX,SYMBL+PATH-1,SYMBL,SYMBL,SYMBL)

		end do
	end do

end program 

!行列の掛け算を行う

subroutine CmplxMultiply(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
	integer A_ROW,A_COL,B_ROW,B_COL
	complex A(A_ROW,A_COL), B(B_ROW,B_COL), C(A_ROW,B_COL)

	if(A_COL.ne.B_ROW) then
		print *, "can't calculate (Multiply)"
		stop
	end if

	do i=1, A_ROW
		do j=1, B_COL
			do k=1, A_COL
				C(i,j) = C(i,j)+A(i,k)*B(k,j)
			end do
		end do
	end do

end subroutine

!行列の和を取る

subroutine CmplxAdd(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
	integer A_ROW,A_COL,B_ROW,B_COL
	complex A(A_ROW,A_COL),B(B_ROW,B_COL),C(A_ROW,A_ROW)

	if((A_ROW.ne.B_ROW).or.(A_COL.ne.B_COL)) then
		print *, "can't calculate (Add)"
		stop
	end if

	do i=1, A_ROW
		do j=1, A_COL
			C(i,j) = A(i,j) + B(i,j)
		end do
	end do 

end subroutine

!随伴行列を取る

subroutine CmplxAdjoint(A,AH,A_ROW,A_COL)
	integer A_ROW,A_COL
	complex A(A_ROW,A_COL),AH(A_COL,A_ROW)

	do i=1, A_ROW
		do j=1, A_COL
			AH(j,i) = conjg(A(i,j))
		end do
	end do

end subroutine

