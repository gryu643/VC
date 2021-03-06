program colReal_noiseNone
	implicit none

	integer,parameter :: H_ROW=1024,H_COL=1024,X_ROW=1024,X_COL=1024,H_PATH=64
	integer,parameter :: NLOOP=200
	integer i,j,k,l,m
	integer SYMBL,PATH
	character(10) TMP
	complex(kind(0d0)) :: Z(H_ROW,H_COL)=(0.0,0.0)
	complex(kind(0d0)) :: X(X_ROW,X_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: Xpre(X_ROW,X_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: H(H_ROW+H_PATH,H_COL)=(0.0,0.0) !(SYMBL+PATH-1,SYMBL)
	complex(kind(0d0)) :: HE(H_ROW+H_PATH*2,H_COL)=(0.0,0.0) !(SYMBL+2*PATH-1,SYMBL+PATH-1)
	complex(kind(0d0)) :: HH(H_ROW,H_COL+H_PATH)=(0.0,0.0) !(SYMBL,SYMBL+PATH-1)
	complex(kind(0d0)) :: HHH(H_ROW,H_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: HX(H_ROW+H_PATH,X_COL)=(0.0,0.0) !(SYMBL+PATH-1,SYMBL)
	complex(kind(0d0)) :: HXarI(H_ROW+H_PATH,X_COL)=(0.0,0.0) !(SYMBL+PATH-1,SYMBL)
	complex(kind(0d0)) :: HHHXbrJ(H_ROW+H_PATH*2,X_COL)=(0.0,0.0) !(SYMBL+2*PATH-1,SYMBL)
	complex(kind(0d0)) :: HHHX(X_ROW,X_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: LAMBDA(X_COL,1)=(0.0,0.0) !(SYMBL,1)
	complex(kind(0d0)) :: Xn(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex(kind(0d0)) :: U(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex(kind(0d0)) :: UH(1,X_ROW)=(0.0,0.0) !(1,SYMBL)
	complex(kind(0d0)) :: LUUH(X_ROW,X_ROW)=(0.0,0.0) !(SYMBL)
	complex(kind(0d0)) :: LUUH_SET(X_ROW,X_ROW)=(0.0,0.0) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: LU(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex(kind(0d0)) :: LUUHXn(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex(kind(0d0)) :: SUB_PART(X_ROW,X_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: arSUB(X_ROW,X_COL)=(0.0,0.0) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: NORM(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex(kind(0d0)) :: UHN(1,1)=(0.0)
	complex(kind(0d0)) :: N(X_ROW,1)=(0.0,0.0) !(SYMBL,1)
	complex(kind(0d0)) :: NAISEKI(1,1)=(0.0,0.0) !(1,1)
	complex(kind(0d0)) :: LAMBDA_MATRIX(X_ROW,X_COL) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: XH(X_ROW,X_COL) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: XLM(X_ROW,X_COL) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: XLMXH(X_ROW,X_COL) !(SYMBL,SYMBL)
	complex(kind(0d0)) :: S(X_ROW,X_COL) !(SYMBL,SYMBL)
	double precision :: LAMBDA_TMP=0.0
	double precision :: TMP1(X_ROW)=0.0 !(SYMBL,1)
	double precision :: NAISEKI_TMP=0.0
	double precision :: sC2=0.0
	double precision :: M1=0.0
	double precision :: M2=0.0
	double precision :: Q(NLOOP,1) !(NLOOP,1)
	double precision :: AVGOTH(NLOOP,1) !(NLOOP,1)

	!シンボル数、パス数を読み込む
	read(5,*) TMP,SYMBL
	read(5,*) TMP,PATH

	!伝搬路行列Hの設定
	do j=0, PATH-1
		do i=1, SYMBL
			H(i+j,i) = cmplx(0.1+0.1*j, 0.2+0.1*j, kind(0d0))
		end do
	end do

	!伝搬路行列Hを拡張したHEを設定
	do j=0, PATH-1
		do i=1, SYMBL+PATH-1
			HE(i+j,i) = cmplx(0.1+0.1*j, 0.2+0.1*j, kind(0d0))
		end do
	end do

	!行列Hの随伴行列HHの設定
	call CAdjoint(H,HH,SYMBL+PATH-1,SYMBL)

	!合成チャネル行列HHHの設定
	call CMultiply(HH,H,HHH,SYMBL,SYMBL+PATH-1,SYMBL+PATH-1,SYMBL)

	do l=1, NLOOP
		!任意伝送ベクトルの設定
		do i=1, SYMBL
			do j=1, SYMBL
				X(i,j) = cmplx(1.0, 0.0, kind(0d0))
			end do
		end do

		do m=1, l
			!次ループでの固有値算出のため、Xを退避
			call CSubstitute(Xpre,X,SYMBL,SYMBL)

			!伝搬路H通過
			call CMultiply(H,X,HX,SYMBL+PATH-1,SYMBL,SYMBL,SYMBL)

			!処理I
			call ProcI(HX,HXarI,SYMBL+PATH-1,SYMBL)	

			!HE通過
			call CMultiply(HE,HXarI,HHHXbrJ,SYMBL+2*(PATH-1),SYMBL+PATH-1,SYMBL+PATH-1,SYMBL)

			!処理J
			call ProcJ(HHHXbrJ,HHHX,SYMBL+2*(PATH-1),SYMBL,PATH)

			!列ベクトル群の固有値をそれぞれ算出
			do i=1, SYMBL
				!列ベクトル内の各行ごとに固有値を算出
				do j=1, SYMBL
					TMP1(j) = abs(real(HHHX(j,i))/real(Xpre(j,i)))
				end do

				!上で出した固有値のうち最小のものをLAMBDA_TMPに格納
				LAMBDA_TMP = TMP1(1)
				do j=2, SYMBL
					if (LAMBDA_TMP.gt.TMP1(j)) then
						LAMBDA_TMP=TMP1(j)
					end if
				end do

				LAMBDA(i,1) = cmplx(LAMBDA_TMP,0.0,kind(0d0))
			end do
			!この時点で配列LAMBDAの各行に固有値が入っている。


			!減算部分の算出
			call CSubstitute(LUUH_SET,Z,SYMBL,SYMBL)
			do i=2, SYMBL
				!収束する固有ベクトル(SYMBL,1)
				do k=1, SYMBL
					Xn(k,1) = Xpre(k,i)
				end do

				!減算する固有ベクトル(SYMBL,1)
				do k=1, SYMBL
					U(k,1) = Xpre(k,i-1)
				end do

				!固有ベクトルの随伴行列(1,SYMBL)
				call CAdjoint(U,UH,SYMBL,1)

				!λ*U(SYMBL,1)
				do k=1, SYMBL
					LU(k,1) = LAMBDA(i-1,1)*U(k,1)
				end do

				!LU*UH(SYMBL,SYMBL)
				call CMultiply(LU,UH,LUUH,SYMBL,1,1,SYMBL)

				!λUUHの集合を格納
				do j=1, SYMBL
					do k=1, SYMBL
						LUUH_SET(j,k) = LUUH_SET(j,k) + LUUH(j,k)
					end do
				end do

				!LUUH_SET*Xn(SYMBL,1) iの次ループで足し合わせる
				call CMultiply(LUUH_SET,Xn,LUUHXn,SYMBL,SYMBL,SYMBL,1)

				!減算部の格納
				do k=1, SYMBL
					SUB_PART(k,i) = LUUHXn(k,1)
				end do

			end do

			!減算
			call CSubtract(HHHX,SUB_PART,arSUB,SYMBL,SYMBL,SYMBL,SYMBL)


			!正規化
			do i=1, SYMBL
				!固有ベクトル群を1列のベクトルに格納
				do j=1, SYMBL
					NORM(j,1) = arSUB(j,i)
				end do

				!正規化
				call CNormalize(NORM,SYMBL,1)

				!正規化したベクトルをXに格納
				do j=1, SYMBL
					X(j,i) = NORM(j,1)
				end do
			end do

		end do

		
		!固有ベクトルか確認(内積=0)
		NAISEKI(1,1) = cmplx(0.0,0.0,kind(0d0))
		do i=1, SYMBL
			do j=i+1, SYMBL
				!固有ベクトル群を１列のベクトルに格納
				do k=1, SYMBL
					U(k,1) = X(k,i)
				end do
				!内積を取る固有ベクトルを格納
				do k=1, SYMBL
					N(k,1) = X(k,j)
				end do

				!随伴行列
				call CAdjoint(U,UH,SYMBL,1)

				!内積の計算
				call CMultiply(UH,N,UHN,1,SYMBL,SYMBL,1)

				!計算した内積を足し合わせる
				call CAdd(NAISEKI,UHN,NAISEKI,1,1,1,1)
			end do
		end do

		!内積の絶対値を出力
		call CAbs(NAISEKI(1,1),NAISEKI_TMP)
		sC2 = SYMBL*(SYMBL-1.0)/2.0
		AVGOTH(l,1) = NAISEKI_TMP / sC2




		!検証
		!スペクトル定理
		do i=1, SYMBL
			LAMBDA_MATRIX(i,i) = LAMBDA(i,1)
		end do
		!スペクトル定理よりシミュレーションで求めた合成チャネル行列を計算
		call CMultiply(X,LAMBDA_MATRIX,XLM,SYMBL,SYMBL,SYMBL,SYMBL)
		call Cadjoint(X,XH,SYMBL,SYMBL)
		call CMultiply(XLM,XH,XLMXH,SYMBL,SYMBL,SYMBL,SYMBL)

		!理論値とシミュレーション結果の差をとる
		call CSubtract(HHH,XLMXH,S,SYMBL,SYMBL,SYMBL,SYMBL)

		!上で計算した差の絶対値の2乗を理論値の絶対値の2乗で正規化する
		M1=0.0
		M2=0.0
		do j=1, SYMBL
			do i=1,SYMBL
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
		complex(kind(0d0)) A(:,:), B(:,:),C(:,:)

		if(A_COL.ne.B_ROW) then
			print *, "can't calculate (Multiply)"
			stop
		end if

		do i=1, A_ROW
			do j=1, B_COL
				C(i,j) = cmplx(0.0,0.0,kind(0d0))
				do k=1, A_COL
					C(i,j) = C(i,j)+A(i,k)*B(k,j)
				end do
			end do
		end do

	end subroutine

	!行列の和を取る

	subroutine CAdd(A,B,C,A_ROW,A_COL,B_ROW,B_COL)
		integer A_ROW,A_COL,B_ROW,B_COL,i,j
		complex(kind(0d0)) A(:,:),B(:,:),C(:,:)

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
		complex(kind(0d0)) A(:,:),B(:,:),C(:,:)

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
		complex(kind(0d0)) A(:,:),AH(:,:)

		do i=1, A_ROW
			do j=1, A_COL
				AH(j,i) = conjg(A(i,j))
			end do
		end do

	end subroutine

	!処理Iを加える

	subroutine ProcI(A,AI,A_ROW,A_COL)
		integer A_ROW,A_COL,i,j
		complex(kind(0d0)) A(:,:),AI(:,:)

		do i=1, A_COL
			do j=1, A_ROW
				AI(A_ROW-j+1,i) = conjg(A(j,i))
			end do
		end do

	end subroutine

	!処理Jを加える

	subroutine ProcJ(A,AJ,A_ROW,A_COL,PATH)
		integer A_ROW,A_COL,PATH,i,j
		complex(kind(0d0)) A(:,:),AJ(:,:)

		do i=1, A_COL
			do j=PATH, A_ROW-(PATH-1)
				AJ(A_ROW-(PATH-1)-j+1,i) = conjg(A(j,i))
			end do
		end do

	end subroutine

	!配列を代入する

	subroutine CSubstitute(A,B,A_ROW,A_COL)
		integer A_ROW,A_COL,i,j
		complex(kind(0d0)) A(:,:),B(:,:)

		do i=1, A_ROW
			do j=1, A_COL
				A(i,j) = B(i,j)
			end do
		end do

	end subroutine

	!正規化する

	subroutine CNormalize(A,A_ROW,A_COL)
		integer A_ROW,A_COL,i
		complex(kind(0d0)) A(:,:)
		double precision TMP

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

	!複素数の絶対値をとる

	subroutine CAbs(A,TMP)
		complex(kind(0d0)) A
		double precision TMP

		!実部と虚部の二乗の和を計算
		TMP = 0.0
		TMP = real(A)**2 + aimag(A)**2
		TMP = sqrt(TMP)
	end subroutine


end program 

