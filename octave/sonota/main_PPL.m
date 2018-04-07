clear
Nsybl = 32;
Npath = 8;
Nlp = 500;

for i=1:Npath
  Ampd(i) = sqrt(1/Npath); %Equal Gain
%  Ampd(i) = sqrt(1/(2^(i-1))); %Exp. atten.
end

for i=1:Npath
  Cpath(i) = complex(randn,randn)/sqrt(2)*Ampd(i);
end

%for i=1:Nsybl
%  for j=1:Npath
%    H(i+j-1,i) = Cpath(j);
%  end
%end
for l=0:Npath-1
  for k = 1:Nsybl
    H(k+l,k) = (0.1+0.1*l) + (0.2+0.1*l)*i;
  end
end
% 
%for i=1:Nsybl+Npath-1
%  for j = 1:Npath
%   HE(i+j-1,i) = Cpath(j);
%  end
%end
for l=0:Npath-1
  for k = 1:Nsybl+Npath-1
    HE(k+l,k) = (0.1+0.1*l) + (0.2+0.1*l)*i;
  end
end
%

for l=1:Nsybl
  for k = 1:Nsybl
    X(l,k) = 1.0 + 0.0*i;
  end
end
tic;
V = PPL (H, HE, X, Nsybl, Npath, Nlp);
toc;
# 固有ベクトルか確認
NAISEKI = 0.0 + 0.0*i;
for l=1:Nsybl
  for k=l+1:Nsybl
    # 固有ベクトル群を１列のベクトルに格納
    U(:,1) = V(:,l);

    # 内積を取る固有ベクトルを格納
    N(:,1) = V(:,k);

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
sC2 = Nsybl*(Nsybl-1.0)/2.0;
AO = NAISEKI_TMP / sC2