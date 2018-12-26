module PPLCombmod
	use CALmod
    use PCONmod2
	use CombSubmod
	implicit none
contains
	subroutine PPLComb(H,HE,X,Eig,Nsybl,Npath,EbN0In,RTNum,UseChNum,ConvStandard,ConvSize,BERStandard,V,Pt,MaxPconLoop)
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
		integer MaxPconLoop

		!declaration
		integer i,j,k,l,m
		complex(kind(0d0)) Xpre(Nsybl,Nsybl)
		complex(kind(0d0)) HX(Nsybl+Npath-1,Nsybl)
		complex(kind(0d0)) HXarI(Nsybl+Npath-1,Nsybl)
		complex(kind(0d0)) HHHXbrJ(Nsybl+2*(Npath-1),Nsybl)
		complex(kind(0d0)) HHHX(Nsybl,Nsybl)
		complex(kind(0d0)) LAMBDA(Nsybl,1)
		complex(kind(0d0)) SUB_PART(Nsybl,Nsybl)
		complex(kind(0d0)) arSUB(Nsybl,Nsybl)
		complex(kind(0d0)) S(Nsybl,Nsybl)
		double precision AVGOTH(ConvSize)
        integer Ksybl(ConvSize)
		double precision BER
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
		Xpre=(0.0,0.0)
		HX=(0.0,0.0)
		HXarI=(0.0,0.0)
		HHHXbrJ=(0.0,0.0)
		HHHX=(0.0,0.0)
		LAMBDA=(0.0,0.0)
		SUB_PART=(0.0,0.0)
		arSUB=(0.0,0.0)
		S=(0.0,0.0)
		Eig=0.0d0
        Ksybl=2
		BER=0.0d0
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

			!固有値算出
			call PPLEV(Xpre,HHHX,LAMBDA,Nsybl)

			!減算部分の算出
			call PPLSubPart(Xpre,LAMBDA,SUB_PART,Nsybl)

			!減算
			call CSubtract(HHHX,SUB_PART,arSUB,Nsybl,Nsybl,Nsybl,Nsybl)

			!正規化
			call PPLNormalize(arSUB,X,Nsybl)
            
			!average othogonality of eiven vector
			call CombOth(X,BERFlag,Ksybl,AVGOTH,AVGOTHN,ConvSize,Nsybl)

			!Ksyblが少ないときの直交度基準を更新
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
						!calculate BER
						call CombIdeal(real(LAMBDA(k,1)),EbN0In(j),Pt_TMP1(1,k),BER)
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

							!calculate BER
							call CombIdeal(real(LAMBDA(i,1)),EbN0In(j),Pt_TMP1(1,i),BER)
							
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
					call CombOut(X,LAMBDA,RTNum,Eig,UseChNum,V,Nsybl,ConvSize,j,l)
					ExitFlag = ExitFlag + 1
				endif
			end do
			
			if(l==MaxPconLoop) then
				call CombOutMax(X,LAMBDA,BERFlag,Nsybl,UseChNum,ConvSize,l,V,Eig,Pt,RTNum)
				ExitFlag = ConvSize
			endif
			!judge PPL exit
			if(ExitFlag==ConvSize) exit
		end do
	end subroutine
end module PPLCombmod
