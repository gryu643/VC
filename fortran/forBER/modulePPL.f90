module modulePPL
	implicit none
contains
	function PPL(H,HE,X,Nsybl,Npath,PPLloop)
		implicit none

		!ifdef
		logical :: PRINT_RSLT=.TRUE.

		!return
		complex(kind(0d0)) :: PPL(Nsybl,Nsybl)=(0.0,0.0)

		!argument
		integer Nsybl,Npath,PPLloop
		complex(kind(0d0)) X(Nsybl,Nsybl)
		complex(kind(0d0)) :: H(Nsybl+Npath-1,Nsybl)
		complex(kind(0d0)) :: HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)

		integer i,j,k,l,m
		complex(kind(0d0)) :: Z(Nsybl,Nsybl)=(0.0,0.0)
		complex(kind(0d0)) :: Xpre(Nsybl,Nsybl)=(0.0,0.0) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: HH(Nsybl,Nsybl+Npath-1)=(0.0,0.0) !(Nsybl,Nsybl+Npath-1)
		complex(kind(0d0)) :: HHH(Nsybl,Nsybl)=(0.0,0.0) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: HX(Nsybl+Npath-1,Nsybl)=(0.0,0.0) !(Nsybl+Npath-1,Nsybl)
		complex(kind(0d0)) :: HXarI(Nsybl+Npath-1,Nsybl)=(0.0,0.0) !(Nsybl+Npath-1,Nsybl)
		complex(kind(0d0)) :: HHHXbrJ(Nsybl+2*(Npath-1),Nsybl)=(0.0,0.0) !(Nsybl+2*Npath-1,Nsybl)
		complex(kind(0d0)) :: HHHX(Nsybl,Nsybl)=(0.0,0.0) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: LAMBDA(Nsybl,1)=(0.0,0.0) !(Nsybl,1)
		complex(kind(0d0)) :: Xn(Nsybl,1)=(0.0,0.0) !(Nsybl,1)
		complex(kind(0d0)) :: U(Nsybl,1)=(0.0,0.0) !(Nsybl,1)
		complex(kind(0d0)) :: UH(1,Nsybl)=(0.0,0.0) !(1,Nsybl)
		complex(kind(0d0)) :: LUUH(Nsybl,Nsybl)=(0.0,0.0) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: LUUH_SET(Nsybl,Nsybl)=(0.0,0.0) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: LU(Nsybl,1)=(0.0,0.0) !(Nsybl,1)
		complex(kind(0d0)) :: LUUHXn(Nsybl,1)=(0.0,0.0) !(Nsybl,1)
		complex(kind(0d0)) :: SUB_PART(Nsybl,Nsybl)=(0.0,0.0) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: arSUB(Nsybl,Nsybl)=(0.0,0.0) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: NORM(Nsybl,1)=(0.0,0.0) !(Nsybl,1)
		complex(kind(0d0)) :: UHN(1,1)=(0.0)
		complex(kind(0d0)) :: N(Nsybl,1)=(0.0,0.0) !(Nsybl,1)
		complex(kind(0d0)) :: NAISEKI(1,1)=(0.0,0.0) !(1,1)
		complex(kind(0d0)) :: LAMBDA_MATRIX(Nsybl,Nsybl) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: XH(Nsybl,Nsybl) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: XLM(Nsybl,Nsybl) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: XLMXH(Nsybl,Nsybl) !(Nsybl,Nsybl)
		complex(kind(0d0)) :: S(Nsybl,Nsybl) !(Nsybl,Nsybl)
		double precision :: LAMBDA_TMP=0.0
		double precision :: TMP1(Nsybl)=0.0 !(Nsybl,1)
		double precision :: NAISEKI_TMP=0.0
		double precision :: sC2=0.0
		double precision :: M1=0.0
		double precision :: M2=0.0
		double precision :: Q(PPLloop,1) !(PPLloop,1)
		double precision :: AVGOTH(PPLloop,1) !(PPLloop,1)
		complex(kind(0d0)) :: D(Nsybl,1)=(0.0,0.0)
		complex(kind(0d0)) :: D1(Nsybl,1)=(0.0,0.0)
		double precision :: D_ABS=0.0
		double precision :: D1_ABS=0.0

		!行列Hの随伴行列HHの設定
		call CAdjoint(H,HH,Nsybl+Npath-1,Nsybl)

		!合成チャネル行列HHHの設定
		call CMultiply(HH,H,HHH,Nsybl,Nsybl+Npath-1,Nsybl+Npath-1,Nsybl)

		do l=1, PPLloop
			!次ループでの固有値算出のため、Xを退避
			call CSubstitute(Xpre,X,Nsybl,Nsybl)

			!伝搬路H通過
			call CMultiply(H,X,HX,Nsybl+Npath-1,Nsybl,Nsybl,Nsybl)

			!処理I
			call ProcI(HX,HXarI,Nsybl+Npath-1,Nsybl)

			!HE通過
			call CMultiply(HE,HXarI,HHHXbrJ,Nsybl+2*(Npath-1),Nsybl+Npath-1,Nsybl+Npath-1,Nsybl)

			!処理J
			call ProcJ(HHHXbrJ,HHHX,Nsybl+2*(Npath-1),Nsybl,Npath)

			!列ベクトル群の固有値をそれぞれ算出
			do i=1, Nsybl
				!列ベクトル内の各行ごとに固有値を算出
				do j=1,Nsybl
					D(j,1) = Xpre(j,i)
					D1(j,1) = HHHX(j,i)
				end do

				call CAbs(D,D_ABS,Nsybl,1)
				call CAbs(D1,D1_ABS,Nsybl,1)

				LAMBDA_TMP = D1_ABS / D_ABS
				LAMBDA(i,1) = cmplx(LAMBDA_TMP,0.0, kind(0d0))
			end do
			!この時点で配列LAMBDAの各行に固有値が入っている。


			!減算部分の算出
			call CSubstitute(LUUH_SET,Z,Nsybl,Nsybl)
			do i=2, Nsybl
				!収束する固有ベクトル(Nsybl,1)
				do k=1, Nsybl
					Xn(k,1) = Xpre(k,i)
				end do

				!減算する固有ベクトル(Nsybl,1)
				do k=1, Nsybl
					U(k,1) = Xpre(k,i-1)
				end do

				!固有ベクトルの随伴行列(1,Nsybl)
				call CAdjoint(U,UH,Nsybl,1)

				!λ*U(Nsybl,1)
				do k=1, Nsybl
					LU(k,1) = LAMBDA(i-1,1)*U(k,1)
				end do

				!LU*UH(Nsybl,Nsybl)
				call CMultiply(LU,UH,LUUH,Nsybl,1,1,Nsybl)

				!λUUHの集合を格納
				do j=1, Nsybl
					do k=1, Nsybl
						LUUH_SET(j,k) = LUUH_SET(j,k) + LUUH(j,k)
					end do
				end do

				!LUUH_SET*Xn(Nsybl,1) iの次ループで足し合わせる
				call CMultiply(LUUH_SET,Xn,LUUHXn,Nsybl,Nsybl,Nsybl,1)

				!減算部の格納
				do k=1, Nsybl
					SUB_PART(k,i) = LUUHXn(k,1)
				end do

			end do

			!減算
			call CSubtract(HHHX,SUB_PART,arSUB,Nsybl,Nsybl,Nsybl,Nsybl)


			!正規化
			do i=1, Nsybl
				!固有ベクトル群を1列のベクトルに格納
				do j=1, Nsybl
					NORM(j,1) = arSUB(j,i)
				end do

				!正規化
				call CNormalize(NORM,Nsybl,1)

				!正規化したベクトルをXに格納
				do j=1, Nsybl
					X(j,i) = NORM(j,1)
				end do
			end do
		end do

		! return
		call CSubstitute(PPL,X,Nsybl,Nsybl)
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
					C(i,j) = cmplx(0.0,0.0, kind(0d0))
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

		subroutine ProcJ(A,AJ,A_ROW,A_COL,Npath)
			integer A_ROW,A_COL,Npath,i,j
			complex(kind(0d0)) A(:,:),AJ(:,:)

			do i=1, A_COL
				do j=Npath, A_ROW-(Npath-1)
					AJ(A_ROW-(Npath-1)-j+1,i) = conjg(A(j,i))
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

		!print

		subroutine print(A)
			complex(kind(0d0)) A(:,:)
			integer i

			do i=1, Nsybl
				print *, l, A(i,1)
			end do
		end subroutine

		!複素数の絶対値をとる

		subroutine CAbs(A,TMP,A_ROW,A_COL)
			integer i,A_ROW,A_COL
			complex(kind(0d0)) A(:,:)
			double precision TMP

			TMP = 0.0
			do i=1, A_ROW
				!実部と虚部の二乗の和を計算
				TMP = TMP + real(A(i,A_COL))**2 + aimag(A(i,A_COL))**2
			end do
			TMP = sqrt(TMP)
		end subroutine

	end function
end module modulePPL