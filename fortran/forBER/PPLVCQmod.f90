module PPLVCQmod
	use CALmod
	implicit none
contains
	subroutine PPLVCQ(H,HE,X,Eig,Nsybl,Npath,RTNum,ConvStandard)
		implicit none

		!argument
		integer Nsybl,Npath
        double precision EbN0In
		complex(kind(0d0)) X(Nsybl,Nsybl)
		complex(kind(0d0)) H(Nsybl+Npath-1,Nsybl)
		complex(kind(0d0)) HE(Nsybl+2*(Npath-1),Nsybl+Npath-1)
		double precision Eig(1,Nsybl)
		double precision ConvStandard
		integer RTNum

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
		double precision AVGOTH
		double precision BER
		complex(kind(0d0)) V1(Nsybl),V2(Nsybl)

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
		BER=0.0d0

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
			do i=1, Nsybl-1
				do j=i+1, Nsybl
					!固有ベクトル群を１列のベクトルに格納
					do k=1, Nsybl
						V1(k) = X(k,i)
					end do
					!内積を取る固有ベクトルを格納
					do k=1, Nsybl
						V2(k) = X(k,j)
					end do

					NAISEKI(1,1) = NAISEKI(1,1) + abs(dot_product(V1,V2))
				end do
			end do

			call CAbs(NAISEKI,NAISEKI_TMP,1,1)
			AVGOTH = NAISEKI_TMP

			if(AVGOTH>ConvStandard) then
				cycle
			else
                do i=1, Nsybl
                    Eig(1,i) = real(LAMBDA(i,1))
                end do
				RTNum = l
                exit
			endif
		end do
	end subroutine
end module PPLVCQmod
