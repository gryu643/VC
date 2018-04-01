function V = PPL (H, HE, X, Nsybl, Npath, lp)
  # Hの随伴行列HHの設定
  HH = H';

  # 合成チャネル豪列HHHの設定
  HHH = HH * H;

  for m=1:lp
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
    for l=Npath:Nsybl+Npath-1
      HHHX(l-(Npath-1),:) = HHHXbrJ(l,:);
    end

    # 列ベクトル群の固有値をそれぞれ算出
    for l=1:Nsybl
      # 列ベクトル内の各行ごとに固有値を算出
      D(:,1) = Xpre(:,l);
      D1(:,1) = HHHX(:,l);

      D_NORM = norm(D);
      D1_NORM = norm(D1);

      LAMBDA_TMP = D1_NORM / D_NORM;
      LAMBDA(l) = LAMBDA_TMP + 0.0*i;
    end

    # 減算部分の導出
    for l=1:Nsybl
      for k=1:Nsybl
        LUUH_SET(k,l) = 0.0 + 0.0*i;
      end
    end

    for l=2:Nsybl
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
    for l=1:Nsybl
      # 固有ベクトル群を1列に格納
      NORM(:,1) = arSUB(:,l);

      # 正規化
      NORM(:,1) = NORM(:,1) / norm(NORM);

      # 正規化したベクトルをXに格納
      X(:,l) = NORM(:,1);
    end
  end

  # 結果の出力
  V = X; 
endfunction