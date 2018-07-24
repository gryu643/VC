function V = PPL (H, HE, X, Nsybl, Npath, lp)
  
  % HH
  HH = H';

  % HHH
  HHH = HH * H;
  for m=1:lp
    % save X
    Xpre = X;

    % pass H
    HX = H * X;

    % procI
    HXarI = conj(HX);
    HXarI = flipud(HXarI);

    % pass HE
    HHHXbrJ = HE * HXarI;

    % procJ
    HHHX = zeros(Nsybl);
    
    HHHXbrJ = conj(HHHXbrJ);
    HHHXbrJ = flipud(HHHXbrJ);
    l = Npath:Nsybl+Npath-1;
    HHHX(l-(Npath-1),:) = HHHXbrJ(l,:);

    
    % eigenvalueD_SUM = zeros(Nsybl,1);
    D_SUM = zeros(Nsybl,1);
    D1_SUM = zeros(Nsybl,1);
    D_NORM = zeros(Nsybl,1);
    D1_NORM = zeros(Nsybl,1);
    LAMBDA_TMP = zeros(Nsybl,1);
    LAMBDA = zeros(Nsybl,1);
    
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
    LAMBDA(l) = LAMBDA_TMP(l);

    % subtraction part
    LUUH_SET = zeros(Nsybl,Nsybl);

    for l=2:Nsybl
      SUB_PART = zeros(Nsybl);
        
      Xn(:,1) = Xpre(:,l);

      U(:,1) = Xpre(:,l-1);

      UH = U';

      LU = LAMBDA(l-1)*U;

      LUUH = LU * UH;
      
      LUUH_SET = LUUH_SET + LUUH;

      LUUHXn = LUUH_SET * Xn;
      
      SUB_PART(:,l) = LUUHXn(:,1);
    end
    arSUB = HHHX - SUB_PART;

%    % subtraction part acceleration ver
%    LUUH = zeros(Nsybl,Nsybl);
%    l=2:Nsybl;
%    UH = Xpre';
%    LU = LAMBDA(l-1).*Xpre(:,l-1);
%    LUUH(:,:,l-1) = LU(:,l-1)*UH(l-1,:);
%    LUUH(:,:,l) = LUUH(:,:,l) + LUUH(:,:,l-1);
%    LUUHXn(l) = LUUH(:,:,l)*Xpre(:,l);
%    arSUB = HHHX - LUUHXn;


    % normalization
    brNORM = zeros(Nsybl);
    NORM = zeros(Nsybl,1);
    
    l=1:Nsybl;
    brNORM(:,l) = arSUB(:,l);
    brNORM = real(brNORM).^2+imag(brNORM).^2;

    brNORM(l) = sum(brNORM,1);

    NORM(l) = sqrt(brNORM(l));

    arSUB(:,l) = arSUB(:,l)./NORM(l);

    X(:,l) = arSUB(:,l);
  end

  % output
  V = X; 
end