program colMat_noiseExist
	implicit none

	integer,parameter :: NLOOP=10
	integer,parameter :: TRIALS=100
	integer i,j,k,l,m,w,ll,seedsize
	integer,allocatable :: seed1(:)
	integer,allocatable :: seed2(:)
	integer SYMBL,PATH
	character(10) TMP
	complex(kind(0d0)) :: Z(100,100)=(0.0,0.0) !(100,100)
	complex(kind(0d0)),allocatable :: X(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: Xpre(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: H(:,:) !(SYMBL+PATH-1,SYMBL)
	complex(kind(0d0)),allocatable :: HE(:,:) !(SYMBL+2*PATH-1,SYMBL+PATH-1)
	complex(kind(0d0)),allocatable :: HH(:,:) !(SYMBL,SYMBL+PATH-1)
	complex(kind(0d0)),allocatable :: HHH(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: HX(:,:) !(SYMBL+PATH-1,SYMBL)
	complex(kind(0d0)),allocatable :: HXarI(:,:) !(SYMBL+PATH-1,SYMBL)
	complex(kind(0d0)),allocatable :: HHHXbrJ(:,:) !(SYMBL+2*PATH-1,SYMBL)
	complex(kind(0d0)),allocatable :: HHHX(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: LAMBDA(:,:) !(SYMBL,1)
	complex(kind(0d0)),allocatable :: Xn(:,:) !(SYMBL,1)
	complex(kind(0d0)),allocatable :: U(:,:) !(SYMBL,1)
	complex(kind(0d0)),allocatable :: UH(:,:) !(1,SYMBL)
	complex(kind(0d0)),allocatable :: LUUH(:,:) !(SYMBL)
	complex(kind(0d0)),allocatable :: LUUH_SET(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: LU(:,:) !(SYMBL,1)
	complex(kind(0d0)),allocatable :: LUUHXn(:,:) !(SYMBL,1)
	complex(kind(0d0)),allocatable :: SUB_PART(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: arSUB(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: NORM(:,:) !(SYMBL,1)
	complex(kind(0d0)) :: UHN(1,1)=(0.0)
	complex(kind(0d0)),allocatable :: N(:,:) !(SYMBL,1)
	complex(kind(0d0)) :: NAISEKI(1,1)=(0.0,0.0) !(1,1)
	complex(kind(0d0)),allocatable :: LAMBDA_MATRIX(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: XH(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: XLM(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: XLMXH(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: S(:,:) !(SYMBL,SYMBL)
	complex(kind(0d0)),allocatable :: Xbefore(:,:) !(SYMBL,SYMBL)
	double precision :: LAMBDA_TMP=0.0
	double precision,allocatable :: TMP1(:) !(SYMBL)
	double precision :: NAISEKI_TMP=0.0
	double precision :: sC2=0.0
	double precision :: M1=0.0
	double precision :: M2=0.0
	double precision :: Q(NLOOP,TRIALS)=0.0 !(NLOOP,TRIALS)
	double precision :: AVGOTH(NLOOP,TRIALS)=0.0 !(NLOOP,TRIALS)
	complex(kind(0d0)),allocatable :: D(:,:) !(SYMBL,1)
	complex(kind(0d0)),allocatable :: D1(:,:) !(SYMBL,1)
	double precision :: D_ABS=0.0
	double precision :: D1_ABS=0.0
	double precision :: EbN0=dble(10.0)
	double precision :: ALPHA=0.0
	complex(kind(0d0)),allocatable :: NOISE2(:,:) !(SYMBL+2*(PATH-1),SYMBL)
	complex(kind(0d0)),allocatable :: NOISE5(:,:) !(SYMBL+PATH-1,SYMBL)
	complex(kind(0d0)),allocatable :: NOISE6(:,:) !(SYMBL+PATH-1,SYMBL)
	complex(kind(0d0)),allocatable :: NOISE7(:,:) !(SYMBL+2*(PATH-1),SYMBL)
	complex(kind(0d0)),allocatable :: NOISE8(:,:) !(SYMBL+2*(PATH-1),SYMBL)
	double precision :: N0=(0.0,0.0)
	complex(kind(0d0)),allocatable :: Bpre(:,:) !(SYMBL+PATH-1,SYMBL)
	complex(kind(0d0)),allocatable :: Bnew(:,:) !(SYMBL+PATH-1,SYMBL)
	double precision :: Qave(NLOOP,1)=0.0
	double precision :: AVGOTHave(NLOOP,1)=0.0
	double precision :: EB=dble(0.0)
	double precision :: EBS=dble(0.0)
	double precision :: Z2=dble(0.0)

	!シンボル数、パス数を読み込む
	read(5,*) TMP,SYMBL
	read(5,*) TMP,PATH
	read(5,*) TMP,ALPHA

	!配列の動的割り付け
	allocate(TMP1(SYMBL))
	allocate(D(SYMBL,1))
	allocate(D1(SYMBL,1))
	allocate(X(SYMBL,SYMBL))
	allocate(Xpre(SYMBL,SYMBL))
	allocate(H(SYMBL+PATH-1,SYMBL))
	allocate(HE(SYMBL+2*(PATH-1),SYMBL+PATH-1))
	allocate(HH(SYMBL,SYMBL+PATH-1))
	allocate(HHH(SYMBL,SYMBL))
	allocate(HX(SYMBL+PATH-1,SYMBL))
	allocate(HXarI(SYMBL+PATH-1,SYMBL))
	allocate(HHHXbrJ(SYMBL+2*(PATH-1),SYMBL))
	allocate(HHHX(SYMBL,SYMBL))
	allocate(LAMBDA(SYMBL,1))
	allocate(Xn(SYMBL,1))
	allocate(U(SYMBL,1))
	allocate(UH(1,SYMBL))
	allocate(LUUH(SYMBL,SYMBL))
	allocate(LUUH_SET(SYMBL,SYMBL))
	allocate(LU(SYMBL,1))
	allocate(LUUHXn(SYMBL,1))
	allocate(SUB_PART(SYMBL,SYMBL))
	allocate(arSUB(SYMBL,SYMBL))
	allocate(NORM(SYMBL,1))
	allocate(N(SYMBL,1))
	allocate(LAMBDA_MATRIX(SYMBL,SYMBL))
	allocate(XH(SYMBL,SYMBL))
	allocate(XLM(SYMBL,SYMBL))
 	allocate(XLMXH(SYMBL,SYMBL))
	allocate(S(SYMBL,SYMBL))
	allocate(Xbefore(SYMBL,SYMBL))
	allocate(NOISE2(SYMBL+2*(PATH-1),SYMBL))
	allocate(NOISE5(SYMBL+PATH-1,SYMBL))
	allocate(NOISE6(SYMBL+PATH-1,SYMBL))
	allocate(NOISE7(SYMBL+2*(PATH-1),SYMBL))
	allocate(NOISE8(SYMBL+2*(PATH-1),SYMBL))
	allocate(Bpre(SYMBL+PATH-1,SYMBL))
	allocate(Bnew(SYMBL+PATH-1,SYMBL))

	!配列初期化
	TMP1=dble(0.0)
	D=cmplx(0.0,0.0,kind(0d0))
	D1=cmplx(0.0,0.0,kind(0d0))
	X=cmplx(0.0,0.0,kind(0d0))
	Xpre=cmplx(0.0,0.0,kind(0d0))
	H=cmplx(0.0,0.0,kind(0d0))
	HE=cmplx(0.0,0.0,kind(0d0))
	HH=cmplx(0.0,0.0,kind(0d0))
	HHH=cmplx(0.0,0.0,kind(0d0))
	HX=cmplx(0.0,0.0,kind(0d0))
	HXarI=cmplx(0.0,0.0,kind(0d0))
	HHHXbrJ=cmplx(0.0,0.0,kind(0d0))
	HHHX=cmplx(0.0,0.0,kind(0d0))
	LAMBDA=cmplx(0.0,0.0,kind(0d0))
	Xn=cmplx(0.0,0.0,kind(0d0))
	U=cmplx(0.0,0.0,kind(0d0))
	UH=cmplx(0.0,0.0,kind(0d0))
	LUUH=cmplx(0.0,0.0,kind(0d0))
	LUUH_SET=cmplx(0.0,0.0,kind(0d0))
	LU=cmplx(0.0,0.0,kind(0d0))
	LUUHXn=cmplx(0.0,0.0,kind(0d0))
	SUB_PART=cmplx(0.0,0.0,kind(0d0))
	arSUB=cmplx(0.0,0.0,kind(0d0))
	NORM=cmplx(0.0,0.0,kind(0d0))
	N=cmplx(0.0,0.0,kind(0d0))
	LAMBDA_MATRIX=cmplx(0.0,0.0,kind(0d0))
	XH=cmplx(0.0,0.0,kind(0d0))
	XLM=cmplx(0.0,0.0,kind(0d0))
 	XLMXH=cmplx(0.0,0.0,kind(0d0))
	S=cmplx(0.0,0.0,kind(0d0))
	Xbefore=cmplx(0.0,0.0,kind(0d0))
	Bpre=cmplx(0.0,0.0,kind(0d0))
	Bnew=cmplx(0.0,0.0,kind(0d0))

	!seedサイズの取得・動的配列の再宣言
	call random_seed(size=seedsize)
	allocate(seed1(seedsize))
	allocate(seed2(seedsize))

	!初期シード設定（Normal関数の中の一様乱数が互いに独立になるため）
	do i=1, seedsize
		seed1(i) = 1
	end do

	do i=1, seedsize
		seed2(i) = 1000000
	end do

	do w=1, TRIALS
		!雑音の設定
		call Setnoise(NOISE2,SYMBL+2*(PATH-1),SYMBL)
		call Setnoise(NOISE5,SYMBL+PATH-1,SYMBL)
		call Setnoise(NOISE6,SYMBL+PATH-1,SYMBL)
		call Setnoise(NOISE7,SYMBL+2*(PATH-1),SYMBL)
		call Setnoise(NOISE8,SYMBL+2*(PATH-1),SYMBL)

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

		ll=1
		do l=1, NLOOP
			if(l.eq.1)then
				ll = 1
			else
				ll = 50*(l-1)
			end if

			!任意伝送ベクトルの設定
			do i=1, SYMBL
				do j=1, SYMBL
					X(i,j) = cmplx(1.0, 0.0, kind(0d0))
				end do
			end do

			do m=1, ll
				!忘却係数適用して固有符号群更新
				do j=1, SYMBL
					do i=1, SYMBL
						X(i,j) = ALPHA*Xbefore(i,j)+(1-ALPHA)*X(i,j)
					end do
				end do

				!次ループでの忘却係数適用のため固有符号群を保持
				call CSubstitute(Xbefore,X,SYMBL,SYMBL)

				!次ループでの固有値算出のため、Xを退避
				call CSubstitute(Xpre,X,SYMBL,SYMBL)

				!増幅値計算
				N0 = dble(0.0)
				do j=1, SYMBL
					do i=1,SYMBL+PATH-1
						N0 = N0 + real(NOISE2(i,j))**2
					end do
				end do
				N0 = N0 / SYMBL

				EB = N0*(10**(EbN0/10))

				EBS = dble(0.0)
				do i=1, SYMBL
					EBS = EBS + real(Xpre(i,1))**2
				end do

				Z2 = sqrt(EB / EBS)

				!雑音の影響を抑えるため増幅
				X = Z2 * X

				!伝搬路H通過
				call CMultiply(H,X,HX,SYMBL+PATH-1,SYMBL,SYMBL,SYMBL)

				!雑音加算
				call Setnoise(NOISE6,SYMBL+PATH-1,SYMBL)
				call CAdd(HX,NOISE6,HX,SYMBL+PATH-1,SYMBL,SYMBL+PATH-1,SYMBL)

				!忘却係数適用
				!今回のループのHXをBnewに格納
				call CSubstitute(Bnew,HX,SYMBL+PATH-1,SYMBL)
				if(m.ne.1) then
					do i=1, SYMBL
						do j=1, SYMBL+PATH-1
							HX(j,i) = ALPHA*Bpre(j,i) + (1-ALPHA)*Bnew(j,i)
						end do
					end do
				end if
				!次ループでの忘却係数適用のため今回のループのHXを格納
				call CSubstitute(Bpre,Bnew,SYMBL+PATH-1,SYMBL)

				!処理I
				call ProcI(HX,HXarI,SYMBL+PATH-1,SYMBL)

				!HE通過
				call CMultiply(HE,HXarI,HHHXbrJ,SYMBL+2*(PATH-1),SYMBL+PATH-1,SYMBL+PATH-1,SYMBL)

				!雑音加算
				call Setnoise(NOISE7,SYMBL+2*(PATH-1),SYMBL)
				call CAdd(HHHXbrJ,NOISE7,HHHXbrJ,SYMBL+2*(PATH-1),SYMBL,SYMBL+2*(PATH-1),SYMBL)
				
				!処理J
				call ProcJ(HHHXbrJ,HHHX,SYMBL+2*(PATH-1),SYMBL,PATH)


				!列ベクトル群の固有値をそれぞれ算出
				do i=1, SYMBL
					!列ベクトル内の各行ごとに固有値を算出
					do j=1,SYMBL
						D(j,1) = Xpre(j,i)
						D1(j,1) = HHHX(j,i)
					end do

					call CAbs(D,D_ABS,SYMBL,1)
					call CAbs(D1,D1_ABS,SYMBL,1)

					LAMBDA_TMP = D1_ABS / D_ABS
					LAMBDA(i,1) = cmplx(LAMBDA_TMP,0.0,kind(0d0))
				end do


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
			call CAbs(NAISEKI,NAISEKI_TMP,1,1)
			sC2 = SYMBL*(SYMBL-1.0)/2.0
			AVGOTH(l,w) = NAISEKI_TMP / sC2


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
			Q(l,w) = M1 / M2

		end do
	end do


	!試行回数分(w=TRIALS)の平均値を計算
	do i=1, TRIALS
		do j=1, NLOOP
			Qave(j,1) = Qave(j,1) + Q(j,i)
			AVGOTHave(j,1) = AVGOTHave(j,1) + AVGOTH(j,i)
		end do
	end do

	do i=1, NLOOP
		Qave(i,1) = Qave(i,1)/dble(TRIALS)
		AVGOTHave(i,1) = AVGOTHave(i,1)/dble(TRIALS)
	end do

	!結果の出力
	do i=1, NLOOP
		print *, i, ",", Qave(i,1), ",", AVGOTHave(i,1)
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

	!正規乱数を配列に設定する

	subroutine Setnoise(A,A_ROW,A_COL)
		integer A_ROW,A_COL,i,j
		complex(kind(0d0)) A(:,:),CNV

		do i=1, A_ROW
			CNV=cmplx(Normal(),Normal(),kind(0d0))

			do j=1, A_COL
				A(i,j) = CNV
			end do
		end do

	end subroutine

	!正規乱数を返す
	function Normal()
		double precision :: m=0.0
		double precision :: var=1.0
		double precision :: x=0.0
		double precision :: y=0.0
		double precision :: pi=dacos(dble(-1))
		double precision :: r1=0.0
		double precision :: r2=0.0
		double precision normal
		double precision :: r=0.0

		!一つ目の一様乱数の生成
		call random_seed(put=seed1)
		call random_number(r1)
		call random_seed(get=seed1)

		!2つ目の一様乱数の生成
		call random_seed(put=seed2)
		call random_number(r2)
		call random_seed(get=seed2)

		!ボックス・ミュラー法による正規乱数生成
		x = sqrt(-2.0*dlog(r1)*var)*dcos(2.0*pi*r2)+m !正規乱数１
		y = sqrt(-2.0*dlog(r1)*var)*dsin(2.0*pi*r2)+m !正規乱数２
		r = sqrt(x**2 + y**2) !レイリー分布

		normal = x
	end function

end program 