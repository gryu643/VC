function V = PPL (H, HE, X, Nsybl, Npath, lp)
  
  # Hの随伴行列HHの設定
  HH = H';

  # 合成チャネル行列HHHの設定
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
    l = Npath:Nsybl+Npath-1;
    HHHX(l-(Npath-1),:) = HHHXbrJ(l,:);
%    for l=Npath:Nsybl+Npath-1
%      HHHX(l-(Npath-1),:) = HHHXbrJ(l,:);
%    end

    # 列ベクトル群の固有値をそれぞれ算出
%    l = 1:Nsybl;
%    D = Xpre;
%    D1 = HHHX;
    
%    D = real(D).^2+imag(D).^2;
%    D1 = real(D1).^2+imag(D1).^2;

%    D_SUM = sum(D,1);
%    D1_SUM = sum(D1,1);

%    D_NORM(l) = sqrt(D_SUM(l));
%    D1_NORM(l) = sqrt(D1_SUM(l));
%
%   LAMBDA_TMP(l) = D1_NORM(l) / D_NORM(l);
%    LAMBDA(l) = LAMBDA_TMP(l) + 0.0*i;
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
    LUUH_SET = zeros(Nsybl,Nsybl);
%    for l=1:Nsybl
%      for k=1:Nsybl
%        LUUH_SET(k,l) = 0.0 + 0.0*i;
%      end
%    end

%    l = 2:Nsybl;
%    LUUHXn = zeros(Nsybl,Nsybl-1);
%    %Xn(Nsybl,Nsybl-1)
%    Xn(:,l-1) = Xpre(:,l);

%    %U(Nsybl,Nsybl-1)
%    U(:,l-1) = Xpre(:,l-1);

%    %UH(Nsybl-1,Nsybl)
%    UH = U';

%    %LU(Nsybl,Nsybl-1)
%    LU(:,l-1) = LAMBDA(l-1).*U(:,l-1);

%    %(1,Nsybl-1)
%    k=1:l-1;
%    LUUHXn(:,l-1) += LU(:,k) * UH(k,:) * Xn(:,l-1);

%    SUB_PART(:,l) = LUUHXn(:,l-1);

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