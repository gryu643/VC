module PPLCombmod
	use CALmod
    use PCONmod2
	implicit none
contains
	subroutine PPLComb(H,HE,X,Eig,Nsybl,Npath,EbN0In,RTNum,UseChNum,ConvStandard,ConvSize,BERStandard,V,Pt)
		implicit none

		!argument
		integer ConvSize,Nsybl,Npath,RTNum(ConvSize),UseChNum(ConvSize)
        double precision EbN0In(ConvSize)
		complex(kind(0d0)) X(Nsybl,Nsybl)
		complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
		complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
		double precision Eig(Nsybl,ConvSize)
		double precision ConvStandard(ConvSize),BERStandard
		complex(kind(0d0)) V(ConvSize,Nsybl,Nsybl)
		double precision Pt(ConvSize,Nsybl)

		!declaration
		integer i,j,k,l,m
		complex(kind(0d0)) Z(Nsybl,Nsybl)
		complex(kind(0d0)) Xpre(Nsybl,Nsybl)
		complex(kind(0d0)) HX(Nsybl+Npath-1,Nsybl)
		complex(kind(0d0)) HXarI(Nsybl+Npath-1,Nsybl)
		complex(kind(0d0)) HHHXbrJ(Nsybl+2*(Npath-1),Nsybl)
		complex(kind(0d0)) HHHX(Nsybl,Nsybl)
		complex(kind(0d0)) LAMBDA(Nsybl,1)
		complex(kind(0d0)) Xn(Nsybl,1)
		complex(kind(0d0)) U(Nsybl,1)
		complex(kind(0d0)) UH(1,Nsybl)
		complex(kind(0d0)) LUUH(Nsybl,Nsybl)
		complex(kind(0d0)) LUUH_SET(Nsybl,Nsybl)
		complex(kind(0d0)) LU(Nsybl,1)
		complex(kind(0d0)) LUUHXn(Nsybl,1)
		complex(kind(0d0)) SUB_PART(Nsybl,Nsybl)
		complex(kind(0d0)) arSUB(Nsybl,Nsybl)
		complex(kind(0d0)) NORM(Nsybl,1)
		complex(kind(0d0)) LAMBDA_MATRIX(Nsybl,Nsybl)
		complex(kind(0d0)) XH(Nsybl,Nsybl)
		complex(kind(0d0)) XLM(Nsybl,Nsybl)
		complex(kind(0d0)) XLMXH(Nsybl,Nsybl)
		complex(kind(0d0)) S(Nsybl,Nsybl)
		double precision LAMBDA_TMP
		complex(kind(0d0)) D(Nsybl,1)
		complex(kind(0d0)) D1(Nsybl,1)
		double precision D_ABS
		double precision D1_ABS
		complex(kind(0d0)) NAISEKI(1,1)
		complex(kind(0d0)) N(Nsybl,1)
		complex(kind(0d0)) UHN(1,1)
		double precision NAISEKI_TMP
		double precision sC2
		double precision AVGOTH(ConvSize)
        integer Ksybl(ConvSize)
		double precision BER
		double precision InstantBER
		double precision LambdaEbN0
		complex(kind(0d0)) V1(Nsybl),V2(Nsybl)
		logical BERFlag(ConvSize)
		integer ExitFlag
		double precision NCS(ConvSize)
		double precision AVGOTHN
		double precision Pt_TMP1(1,Nsybl)
		double precision Pt_TMP2(1,Nsybl)
		integer info
		double precision EbN0Pcon
		double precision EigPcon(1,Nsybl)

		!initialize
		Z=(0.0,0.0)
		Xpre=(0.0,0.0)
		HX=(0.0,0.0)
		HXarI=(0.0,0.0)
		HHHXbrJ=(0.0,0.0)
		HHHX=(0.0,0.0)
		LAMBDA=(0.0,0.0)
		Xn=(0.0,0.0)
		U=(0.0,0.0)
		UH=(0.0,0.0)
		LUUH=(0.0,0.0)
		LUUH_SET=(0.0,0.0)
		LU=(0.0,0.0)
		LUUHXn=(0.0,0.0)
		SUB_PART=(0.0,0.0)
		arSUB=(0.0,0.0)
		NORM=(0.0,0.0)
		LAMBDA_MATRIX=(0.0,0.0)
		XH=(0.0,0.0)
		XLM=(0.0,0.0)
		XLMXH=(0.0,0.0)
		S=(0.0,0.0)
		LAMBDA_TMP=0.0
		D=(0.0,0.0)
		D1=(0.0,0.0)
		D_ABS=0.0
		D1_ABS=0.0
		Eig=0.0d0
        Ksybl=2
		BER=0.0d0
		InstantBER=0.0d0
		LambdaEbN0=0.0d0
		ExitFlag=0
		BERFlag=.False.
		info=1
		Pt_TMP1=0.0d0

		l = 0
		do
			l = l + 1
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
				call CAdd(LUUH_SET,LUUH,LUUH_SET,Nsybl,Nsybl,Nsybl,Nsybl)

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
            
			!average othogonality of eiven vector
			NAISEKI(1,1) = cmplx(0.0,0.0,kind(0d0))
			do i=2, Nsybl
				do j=1, i-1
					!固有ベクトル群を１列のベクトルに格納
					do k=1, Nsybl
						V1(k) = X(k,i)
					end do
					!内積を取る固有ベクトルを格納
					do k=1, Nsybl
						V2(k) = X(k,j)
					end do
					NAISEKI(1,1) = NAISEKI(1,1) + abs(dot_product(V1,V2))

					call CAbs(NAISEKI,NAISEKI_TMP,1,1)
					do m=1, ConvSize
						if(BERFlag(m)) cycle

						if(i==Ksybl(m)) then
							AVGOTH(m) = NAISEKI_TMP
						endif
					end do
				end do
			end do
			call CAbs(NAISEKI,NAISEKI_TMP,1,1)
			AVGOTHN = NAISEKI_TMP

			do k=1, ConvSize
				if(Ksybl(k)<=Nsybl/8) then
					NCS(k) = ConvStandard(k)/10.0d0
				else
					NCS(K) = ConvStandard(k)
				endif
				!NCS(K) = ConvStandard(k)
			end do

			do j=1, ConvSize
				if(BERFlag(j)) cycle

				!judge convergence for Nsybl
				if(AVGOTHN<=ConvStandard(j)) then
					!Transmit power control
					Pt_TMP1=0.0d0
					do k=1, Nsybl
						EigPcon(1,k) = real(LAMBDA(k,1))
					end do
					EbN0Pcon = EbN0In(j)
					call Pcontrol2(EigPcon,EbN0Pcon,Pt_TMP1,Nsybl,Nsybl,info)

					!judge ber
					BER=0.0d0
					do k=1, Nsybl
						LambdaEbN0 = real(LAMBDA(k,1))*EbN0In(j)*Pt_TMP1(1,k)
						InstantBER = 1.0d0/2.0d0*erfc(sqrt(LambdaEbN0))
						BER = BER + InstantBER
					end do
					if(BER/dble(Nsybl)<=BERStandard) then
						!output
						BERFlag(j) = .True.
						UseChNum(j) = Nsybl
						do k=1, Nsybl
							Pt(j,k) = Pt_TMP1(1,k)
						end do
					endif
				endif

				!judge convergence for Ksybl
				if(BERFlag(j).eqv..False.) then
					if(AVGOTH(j)>NCS(j)) then
						cycle
					else
						!judge ber
						BER=0.0d0
						do i=1, Ksybl(j)
							!transmit power control
							!save taransmit Power
							do k=1, Nsybl
								Pt(j,k) = Pt_TMP1(1,k)
							end do
							Pt_TMP1=0.0d0
							do k=1, Nsybl
								EigPcon(1,k) = real(LAMBDA(k,1))
							end do
							EbN0Pcon = EbN0In(j)
							call Pcontrol2(EigPcon,EbN0Pcon,Pt_TMP1,i,Nsybl,info)

							LambdaEbN0 = real(LAMBDA(i,1))*EbN0In(j)*Pt_TMP1(1,i)
							InstantBER = 1.0d0/2.0d0*erfc(sqrt(LambdaEbN0))
							BER = BER + InstantBER
							if(BER/dble(i)>1.0d0/2.0d0*erfc(sqrt(real(LAMBDA(i,1))*EbN0In(j)))/dble(i)) then
								print *, Ksybl(j),i,BER/dble(i), 1.0d0/2.0d0*erfc(sqrt(real(LAMBDA(i,1))*EbN0In(j)))/dble(i)
							endif
							
							!skip judgement to BER of Ksybl when Ksybl>2
							if(Ksybl(j)>2.and.i<Ksybl(j)) cycle

							if(BER/dble(i)>BERStandard) then
								!output and exit
								BERFlag(j) = .True.
								UseChNum(j) = i-1
								exit
							else
								if(Ksybl(j)==2.and.i==1) cycle

								Ksybl(j) = Ksybl(j) + 1
								if(Ksybl(j)==Nsybl+1) then
									!output and exit
									BERFlag(j) = .True.
									UseChNum(j) = Nsybl
									do k=1, Nsybl
										Pt(j,k) = Pt_TMP1(1,k)
									end do
									exit
								endif
							endif
						end do
					endif
				endif

				!output result 
				if(BERFlag(j)) then
					RTNum(j) = l
					ExitFlag = ExitFlag + 1
					do m=1, Nsybl
						Eig(m,j) = real(LAMBDA(m,1))
					end do
			
					select case(UseChNum(j))
						case(0)
						case(1:)
							do m=1, UseChNum(j)
								do k=1, Nsybl
									V(j,k,m) = X(k,m)
								end do
							end do
					end select
				endif
			end do
			
			!judge PPL exit
			if(ExitFlag==ConvSize) exit
		end do
	end subroutine
end module PPLCombmod
