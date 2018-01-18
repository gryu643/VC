program seqMat_noiseNone
	implicit none

	integer,parameter :: H_ROW=1024,H_COL=1024,X_ROW=1024,X_COL=1024,H_PATH=64
	integer,parameter :: NLOOP=200
	integer i,j,k,l,m,n
	integer SYMBL,PATH
	character(10) TMP
	complex :: Z(H_ROW,H_COL)=(0.0,0.0)
	complex :: XG(X_ROW,X_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex :: X(X_ROW,1) !(SYMBL,1)
	complex :: Xpre(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex :: H(H_ROW+H_PATH,H_COL)=(0.0,0.0) !(SYMBL+PATH-1,SYMBL)
	complex :: HE(H_ROW+H_PATH*2,H_COL)=(0.0,0.0) !(SYMBL+2*PATH-1,SYMBL+PATH-1)
	complex :: HH(H_ROW,H_COL+H_PATH)=(0.0,0.0) !(SYMBL,SYMBL+PATH-1)
	complex :: HHH(H_ROW,H_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex :: HX(H_ROW+H_PATH,1)=(0.0,0.0) !(SYMBL+PATH-1,1)
	complex :: HXarI(H_ROW+H_PATH,1)=(0.0,0.0) !(SYMBL+PATH-1,1)
	complex :: HHHXbrJ(H_ROW+H_PATH*2,1)=(0.0,0.0) !(SYMBL+2*PATH-1,1)
	complex :: HHHX(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex :: LAMBDA(X_COL,1)=(0.0,0.0) !(SYMBL,1)
	complex :: Xn(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex :: U(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex :: UH(1,X_ROW)=(0.0,0.0) !(1,SYMBL)
	complex :: LUUH(X_ROW,X_ROW)=(0.0,0.0) !(SYMBL,SYMBL)
	complex :: LUUH_SET(X_ROW,X_ROW)=(0.0,0.0) !(SYMBL,SYMBL)
	complex :: LU(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex :: LUUHXn(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex :: arSUB(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex :: NORM(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex :: UHN(1,1)=(0.0)
	complex :: N1(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex :: NAISEKI(1,1)=(0.0,0.0) !(1,1)
	complex :: LAMBDA_MATRIX(X_ROW,X_COL) !(SYMBL,SYMBL)
	complex :: XGH(X_ROW,X_COL) !(SYMBL,SYMBL)
	complex :: XGLM(X_ROW,X_COL) !(SYMBL,SYMBL)
	complex :: XGLMXGH(X_ROW,X_COL) !(SYMBL,SYMBL)
	complex :: S(X_ROW,X_COL) !(SYMBL,SYMBL)
	real :: LAMBDA_TMP=0.0
	real :: TMP1(X_ROW)=0.0 !(SYMBL,1)
	real :: NAISEKI_TMP=0.0
	real :: sC2=0.0
	real :: M1=0.0
	real :: M2=0.0
	real :: Q(NLOOP,1) !(NLOOP,1)
	real :: AVGOTH(NLOOP,1) !(NLOOP,1)
	complex :: D(X_COL,1)=(0.0,0.0)
	complex :: D1(X_COL,1)=(0.0,0.0)
	real :: D_ABS=0.0
	real :: D1_ABS=0.0

	!シンボル数、パス数を読み込む
	read(5,*) TMP,SYMBL
	read(5,*) TMP,PATH

	!伝搬路行列Hの設定
	do j=0, PATH-1
		do i=1, SYMBL
			H(i+j,i) = cmplx(0.1+0.1*j, 0.2+0.1*j)
		end do
	end do

	!伝搬路行列Hを拡張したHEを設定
	do j=0, PATH-1
		do i=1, SYMBL+PATH-1
			HE(i+j,i) = cmplx(0.1+0.1*j, 0.2+0.1*j)
		end do
	end do

	!行列Hの随伴行列HHの設定
	call CAdjoint(H,HH,SYMBL+PATH-1,SYMBL)

	!合成チャネル行列HHHの設定
	call CMultiply(HH,H,HHH,SYMBL,SYMBL+PATH-1,SYMBL+PATH-1,SYMBL)

	!固有符号ひとつ当たりの往復回数
	do l=1, NLOOP
		!任意伝送ベクトルの設定・初期化
		do i=1, SYMBL
			do j=1, SYMBL
				XG(i,j) = cmplx(1.0, 0.0)
			end do
		end do

		!固有値行列を初期化
		call CSubstitute(LAMBDA,Z,SYMBL,1)

		!減算部のLUUH_SETを初期化
		call CSubstitute(LUUH_SET,Z,SYMBL,SYMBL)

		!固有符号を第1固有符号から順番に収束させる
		do m=1, SYMBL
			do i=1, SYMBL
				X(i,1) = XG(i,m)
			end do

			!各固有符号をl回往復
			do n=1, l
				!固有値算出のため、Xを退避
				call CSubstitute(Xpre,X,SYMBL,1)

				!伝搬路H通過
				call CMultiply(H,X,HX,SYMBL+PATH-1,SYMBL,SYMBL,1) !HX(SYMBL+PATH-1,1)

				!処理I
				call ProcI(HX,HXarI,SYMBL+PATH-1,1) !HXarI(SYMBL+PATH-1,1)

				!HE通過
				call CMultiply(HE,HXarI,HHHXbrJ,SYMBL+2*(PATH-1),SYMBL+PATH-1,SYMBL+PATH-1,1) !HHHXbrJ(SYMBL+2*(PATH-1), 1)

				!処理J
				call ProcJ(HHHXbrJ,HHHX,SYMBL+2*(PATH-1),1,PATH) !HHHX(SYMBL,1)


				!減算処理
				if(m.gt.1) then
					!往復前の固有ベクトルをXnに格納
					call CSubstitute(Xn,Xpre,SYMBL,1)

					!LUUH_SET*Xn(SYMBL,1)
					call CMultiply(LUUH_SET,Xn,LUUHXn,SYMBL,SYMBL,SYMBL,1)

					!減算
					call CSubtract(HHHX,LUUHXn,arSUB,SYMBL,SYMBL,SYMBL,SYMBL)
				else if(m.eq.1) then
					call CSubstitute(arSUB,HHHX,SYMBL,1)
				end if


				!正規化
				!固有ベクトル群を1列のベクトルに格納
				call CSubstitute(NORM,arSUB,SYMBL,1)

				!正規化
				call CNormalize(NORM,SYMBL,1)

				!正規化したベクトルをXに格納
				call CSubstitute(X,NORM,SYMBL,1)
			end do


			!固有符号群行列に求めた固有符号を格納
			do i=1, SYMBL
				XG(i,m) = Xpre(i,1)
			end do


			!固有値の導出
			!列ベクトルの固有値を算出　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　
			call CSubstitute(D,Xpre,SYMBL,1)
			call CSubstitute(D1,HHHX,SYMBL,1)

			call CAbs(D,D_ABS,SYMBL,1)
			call CAbs(D1,D1_ABS,SYMBL,1)

			LAMBDA_TMP = D1_ABS / D_ABS
			LAMBDA(m,1) = cmplx(LAMBDA_TMP,0.0)

			!減算部分の算出
			!減算する固有ベクトル U(SYMBL,1)
			call CSubstitute(U,Xpre,SYMBL,1)

			!固有ベクトルの随伴行列 UH(1,SYMBL)
			call CAdjoint(U,UH,SYMBL,1)

			!λ*U(SYMBL,1)
			do k=1, SYMBL
				LU(k,1) = LAMBDA(m,1)*U(k,1)
			end do

			!LU*UH(SYMBL,SYMBL)
			call CMultiply(LU,UH,LUUH,SYMBL,1,1,SYMBL)

			!λUUHの集合を格納
			do j=1, SYMBL
				do k=1, SYMBL
					LUUH_SET(j,k) = LUUH_SET(j,k) + LUUH(j,k)
				end do
			end do
		end do


		!固有ベクトルか確認(内積=0)
		NAISEKI(1,1) = cmplx(0.0,0.0)
		do i=1, SYMBL
			do j=i+1, SYMBL
				!固有ベクトル群を１列のベクトルに格納
				do k=1, SYMBL
					U(k,1) = XG(k,i)
				end do
				!内積を取る固有ベクトルを格納
				do k=1, SYMBL
					N1(k,1) = XG(k,j)
				end do

				!随伴行列
				call CAdjoint(U,UH,SYMBL,1)

				!内積の計算
				call CMultiply(UH,N1,UHN,1,SYMBL,SYMBL,1)

				!計算した内積を足し合わせる
				call CAdd(NAISEKI,UHN,NAISEKI,1,1,1,1)
			end do
		end do


		!平均直交度を出力
		call CAbs(NAISEKI,NAISEKI_TMP,1,1)
		sC2 = SYMBL*(SYMBL-1.0)/2.0
		AVGOTH(l,1) = NAISEKI_TMP / sC2


		!検証
		!スペクトル定理
		do i=1, SYMBL
			LAMBDA_MATRIX(i,i) = LAMBDA(i,1)
		end do
		!スペクトル定理よりシミュレーションで求めた合成チャネル行列を計算
		call CMultiply(XG,LAMBDA_MATRIX,XGLM,SYMBL,SYMBL,SYMBL,SYMBL)
		call Cadjoint(XG,XGH,SYMBL,SYMBL)
		call CMultiply(XGLM,XGH,XGLMXGH,SYMBL,SYMBL,SYMBL,SYMBL)

		!理論値とシミュレーション結果の差をとる
		call CSubtract(HHH,XGLMXGH,S,SYMBL,SYMBL,SYMBL,SYMBL)

		!上で計算した差の絶対値の2乗を理論値の絶対値の2乗で正規化する
		M1=0.0
		M2=0.0
		do i=1, SYMBL
			do j=1,SYMBL
				M1 = M1 + real(S(i,j))**2 + aimag(S(i,j))**2
				M2 = M2 + real(HHH(i,j))**2 + aimag(HHH(i,j))**2
			end do
		end do
		Q(l,1) = M1 / M2

	end do

	!結果の出力
	do i=1, NLOOP
		print *, i, ",", Q(i,1), ",", AVGOTH(i,1)
	end do

contains
	!行列の掛け算を行う

	subroutine CMultiply(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
		integer A_ROW,A_COL,B_ROW,B_COL,i,j,k
		complex A(:,:), B(:,:),C(:,:)

		if(A_COL.ne.B_ROW) then
			print *, "can't calculate (Multiply)"
			stop
		end if

		do i=1, A_ROW
			do j=1, B_COL
				C(i,j) = cmplx(0.0,0.0)
				do k=1, A_COL
					C(i,j) = C(i,j)+A(i,k)*B(k,j)
				end do
			end do
		end do

	end subroutine

	!行列の和を取る

	subroutine CAdd(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
		integer A_ROW,A_COL,B_ROW,B_COL,i,j
		complex A(:,:),B(:,:),C(:,:)

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

	!行列の差を取る

	subroutine CSubtract(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
		integer A_ROW,A_COL,B_ROW,B_COL,i,j
		complex A(:,:),B(:,:),C(:,:)

		if((A_ROW.ne.B_ROW).or.(A_COL.ne.B_COL)) then
			print *, "can't calculate (Subtract)"
			stop
		end if

		do i=1, A_ROW
			do j=1, A_COL
				C(i,j) = A(i,j) - B(i,j)
			end do
		end do 

	end subroutine

	!随伴行列を取る

	subroutine CAdjoint(A,AH,A_ROW,A_COL)
		integer A_ROW,A_COL,i,j
		complex A(:,:),AH(:,:)

		do i=1, A_ROW
			do j=1, A_COL
				AH(j,i) = conjg(A(i,j))
			end do
		end do

	end subroutine

	!処理Iを加える

	subroutine ProcI(A,AI,A_ROW,A_COL)
		integer A_ROW,A_COL,i,j
		complex A(:,:),AI(:,:)

		do i=1, A_COL
			do j=1, A_ROW
				AI(A_ROW-j+1,i) = conjg(A(j,i))
			end do
		end do

	end subroutine

	!処理Jを加える

	subroutine ProcJ(A,AJ,A_ROW,A_COL,PATH)
		integer A_ROW,A_COL,PATH,i,j
		complex A(:,:),AJ(:,:)

		do i=1, A_COL
			do j=PATH, A_ROW-(PATH-1)
				AJ(A_ROW-(PATH-1)-j+1,i) = conjg(A(j,i))
			end do
		end do

	end subroutine

	!配列を代入する

	subroutine CSubstitute(A,B,A_ROW,A_COL)
		integer A_ROW,A_COL,i,j
		complex A(:,:),B(:,:)

		do i=1, A_ROW
			do j=1, A_COL
				A(i,j) = B(i,j)
			end do
		end do

	end subroutine

	!正規化する

	subroutine CNormalize(A,A_ROW,A_COL)
		integer A_ROW,A_COL,i
		complex A(:,:)
		real TMP

		if(A_COL.ne.1) then
			print *, "The number of column isn't one."
			stop
		end if

		!実部と虚部の二乗の和を計算
		TMP = 0.0
		do i=1, A_ROW
			TMP = TMP + real(A(i,1))**2 + aimag(A(i,1))**2
		end do

		TMP = sqrt(TMP)

		!各成分をTMPで割る
		do i=1, A_ROW
			A(i,1) = A(i,1) / TMP
		end do

	end subroutine

	!print

	subroutine print(A)
		complex A(:,:)
		integer i

		do i=1, SYMBL
			print *, l, A(i,2)
		end do
	end subroutine

	!複素数の絶対値をとる

	subroutine CAbs(A,TMP,A_ROW,A_COL)
		integer i,A_ROW,A_COL
		complex A(:,:)
		real TMP

		TMP = 0.0
		do i=1, A_ROW
			!実部と虚部の二乗の和を計算
			TMP = TMP + real(A(i,A_COL))**2 + aimag(A(i,A_COL))**2
		end do
		TMP = sqrt(TMP)
	end subroutine


end program 