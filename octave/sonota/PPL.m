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

    # 列ベクトル群の固有値をそれぞれ算出
    l = 1:Nsybl;
    D = Xpre;
    D1 = HHHX;
    
    D = real(D).^2+imag(D).^2;
    D1 = real(D1).^2+imag(D1).^2;

    D_SUM(l) = sum(D,1);
    D1_SUM(l) = sum(D1,1);

    D_NORM(l) = sqrt(D_SUM(l));
    D1_NORM(l) = sqrt(D1_SUM(l));

    LAMBDA_TMP(l) = D1_NORM(l) ./ D_NORM(l);
    LAMBDA(l) = LAMBDA_TMP(l) + 0.0*i;

    # 減算部分の導出
    LUUH_SET = zeros(Nsybl,Nsybl);

    l = 2:Nsybl;
    LUUH = zeros(Nsybl,Nsybl,Nsybl);
    SUB_PART = zeros(Nsybl,Nsybl);

    Xn(:,l-1) = Xpre(:,l);

    U(:,l-1) = Xpre(:,l-1);

    UH = U';

    LU(:,l-1) = LAMBDA(l-1).*U(:,l-1);

    for j=2:Nsybl
        k = 1:l(j)
        SUB_PART(:,j) = LU(:,k) * UH(k,:) * Xn(:,j-1);
    end
%    for k=2:Nsybl
%        LUUH(:,:,k) = LU(:,k-1) * UH(k-1,:);
%    end 
%
%    for k=2:Nsybl
%        LUUH(:,:,k) = LUUH(:,:,k-1) + LUUH(:,:,k);
%    end
%    for k=2:Nsybl        
%        SUB_PART(:,k) = LUUH(:,:,k) * Xn(:,k-1);
%    end

%    for l=2:Nsybl
%      # 収束する固有ベクトル
%      Xn(:,1) = Xpre(:,l);
%
%      # 減算する固有ベクトル
%      U(:,1) = Xpre(:,l-1);
%
%      # 固有ベクトルの随伴行列
%      UH = U';
%
%      # λ*U
%      LU = LAMBDA(l-1)*U;
%
%      # LU*UH
%      LUUH = LU * UH;
%
%      # LUUHの集合を格納
%      LUUH_SET = LUUH_SET + LUUH;
%
%      # LUUH_SET*Xn
%      LUUHXn = LUUH_SET * Xn;
%      
%      # 減算部の格納
%      SUB_PART(:,l) = LUUHXn(:,1);
%    end
    arSUB = HHHX - SUB_PART;

    l=1:Nsybl;
    brNORM(:,l) = arSUB(:,l);
    brNORM = real(brNORM).^2+imag(brNORM).^2;

    brNORM(l) = sum(brNORM,1);

    NORM(l) = sqrt(brNORM(l));

    arSUB(:,l) = arSUB(:,l) ./ NORM(l);

    X(:,l) = arSUB(:,l);
  end

  # 結果の出力
  V = X; 
endfunction