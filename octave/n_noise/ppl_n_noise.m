clear;
SYMBL =32;
PATH = 8;
NLOOP = 500;

# 伝搬路行列Hの設定
for l=0:PATH-1
	for k = 1:SYMBL
		H(k+l,k) = (0.1+0.1*l) + (0.2+0.1*l)*i;
	end
end

# 伝搬路行列Hを拡張したHEを設定
for l=0:PATH-1
	for k = 1:SYMBL+PATH-1
		HE(k+l,k) = (0.1+0.1*l) + (0.2+0.1*l)*i;
	end
end

# Hの随伴行列HHの設定
HH = H';

# 合成チャネル豪列HHHの設定
HHH = HH * H;

# 任意伝送ベクトルの設定
for l=1:SYMBL
	for k = 1:SYMBL
		X(l,k) = 1.0 + 0.0*i;
	end
end

for m=1:NLOOP
	# 固有値算出のためX退避
	Xpre = X;

	# 伝搬路H通過
	HX = H * X;

	# 処理I
	HXarI = conj(HX);
	HXarI = flipud(HXarI);

	# HE通過
	HHHXbrJ = HE * HXarI;

	# 処理J
	HHHXbrJ = conj(HHHXbrJ);
	HHHXbrJ = flipud(HHHXbrJ);
	for l=PATH:SYMBL+PATH-1
		HHHX(l-(PATH-1),:) = HHHXbrJ(l,:);
	end

	# 列ベクトル群の固有値をそれぞれ算出
	for l=1:SYMBL
		# 列ベクトル内の各行ごとに固有値を算出
		D(:,1) = Xpre(:,l);
		D1(:,1) = HHHX(:,l);

		D_NORM = norm(D);
		D1_NORM = norm(D1);

		LAMBDA_TMP = D1_NORM / D_NORM;
		LAMBDA(l) = LAMBDA_TMP + 0.0*i;
	end

	# 減算部分の導出
	for l=1:SYMBL
		for k=1:SYMBL
			LUUH_SET(k,l) = 0.0 + 0.0*i;
		end
	end

	for l=2:SYMBL
		# 収束する固有ベクトル
		Xn(:,1) = Xpre(:,l);

		# 減算する固有ベクトル
		U(:,1) = Xpre(:,l-1);

		# 固有ベクトルの随伴行列
		UH = U';

		# λ*U
		LU = LAMBDA(l-1)*U;

		# LU*UH
		LUUH = LU * UH;

		# LUUHの集合を格納
		LUUH_SET = LUUH_SET + LUUH;

		# LUUH_SET*Xn
		LUUHXn = LUUH_SET * Xn;
		
		# 減算部の格納
		SUB_PART(:,l) = LUUHXn(:,1);
	end
	arSUB = HHHX - SUB_PART;

	# 正規化
	for l=1:SYMBL
		# 固有ベクトル群を1列に格納
		NORM(:,1) = arSUB(:,l);

		# 正規化
		NORM(:,1) = NORM(:,1) / norm(NORM);

		# 正規化したベクトルをXに格納
		X(:,l) = NORM(:,1);
	end

	# 固有ベクトルか確認
	NAISEKI = 0.0 + 0.0*i;
	for l=1:SYMBL
		for k=l+1:SYMBL
			# 固有ベクトル群を１列のベクトルに格納
			U(:,1) = X(:,l);

			# 内積を取る固有ベクトルを格納
			N(:,1) = X(:,k);

			# 随伴行列
			UH = U';

			# 内積の計算
			UHN = UH * N;

			# 計算した内積を足し合わせる
			NAISEKI = NAISEKI + UHN;
		end
	end

	# 内積の絶対値を出力
	NAISEKI_TMP = norm(NAISEKI);
	sC2 = SYMBL*(SYMBL-1.0)/2.0;
	AO(m,1) = m;
	AO(m,2) = NAISEKI_TMP / sC2;

	# 検証
	# スペクトル定理
	for l=1:SYMBL
		LAMBDA_MATRIX(l,l) = LAMBDA(l);
	end

	# スペクトル定理よりシミュレーションで求めた合成チャネル行列を計算
	XLM = X * LAMBDA_MATRIX;
	XH = X';
	XLMXH = XLM * XH;

	# 理論値とシミュレーション結果の差を取る
	S = HHH - XLMXH;

	# 上で計算した差の絶対値の２乗を理論値の絶対値の２乗で正規化する
	M1=0.0;
	M2=0.0;
	for l=1:SYMBL
		for k=1:SYMBL
			M1 = M1 + real(S(k,l))^2 + imag(S(k,l))^2;
			M2 = M2 + real(HHH(k,l))^2 + imag(HHH(k,l))^2;
		end
	end
	Q(m,1) = m;
	Q(m,2) = M1 / M2;
end

# 結果の出力
dlmwrite('out.csv',[Q(:,1) Q(:,2)],'delimiter',',');