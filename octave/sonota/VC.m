close all;
clear all;
tic;
Nsybl = 32;
Npath = 8;
SEbN0 = -10;
EEbN0 = 40;
Step = 5;

Nloop = 10000;
PPLloop = 500;

fileID = fopen('ber(Npath=8,Nsym=32).txt','w');


%Channel parameter setting
%for i=1:Npath
%  Cpath(i) = complex((0.1+0.1*(i-1)),(0.2 + 0.1*(i-1)));
%end
%
%Channel gain parameter
i=1:Npath;
Ampd(i) = sqrt(1/Npath);
%%for i=1:Npath
%%  Ampd(i) = sqrt(1/Npath); %Equal Gain
%%%  Ampd(i) = sqrt(1/(2^(i-1))); %Exp. atten.
%%end
%
for KEbN0=SEbN0:Step:EEbN0 %Eb/N0 loop
  Psig = 0;
  Pwgn = 0;
  Collect = 0;
  False = 0;
%
  for loop=1:Nloop %Monte calro loop
%
%Channel parameter setting
    i=1:Npath;
    Cpath(i) = complex(randn,randn)/sqrt(2)*Ampd(i);
%%    for i=1:Npath
%%      Cpath(i) = complex(randn,randn)/sqrt(2)*Ampd(i);
%%    end
%
    for i=1:Nsybl
      j = 1:Npath;
      H(i+j-1,i) = Cpath(j);
    end
% 
    for i=1:Nsybl+Npath-1
      j = 1:Npath;
      HE(i+j-1,i) = Cpath(j);
    end
%
    l = 1:Nsybl;
    Xppl(l,l) = 1.0 + 0.0*i;
%%    for l=1:Nsybl
%%      for k = 1:Nsybl
%%        Xppl(l,k) = 1.0 + 0.0*i;
%%      end
%%    end
%
    V = PPL (H, HE, Xppl, Nsybl, Npath, PPLloop);
%    Xppl=V;
    
    HH = ctranspose(H);
    HHH = HH*H;
%    [V,D] = eig(HHH);
%

    S = complex(round(rand(1,Nsybl))*2-1,round(rand(1,Nsybl))*2-1); %[-1,1] Transmit symbol
    TdatI = (real(S)+1)/2;
    TdatQ = (imag(S)+1)/2;
    i=1:Nsybl;
    SU(:,i) = S(i).*V(:,i);
%%    for i=1:Nsybl;
%%      SU(:,i) = S(i)*V(:,i);
%%    end
%
    X = zeros(1,Nsybl);
    for j=1:Nsybl;
      for i=1:Nsybl;
        X(i) = X(i) + SU(i,j); %Transmit vector
      end
    end
%
    Pow = 0;
    for i=1:Nsybl
      Pow = Pow + abs(X(i))^2/Nsybl;
    end
    X /= sqrt(Pow); %Transmit power control
    Y = H*transpose(X); %Channel
%
    Noise = complex(randn(Nsybl+Npath-1,1),randn(Nsybl+Npath-1,1));
    Noise *= sqrt(1/10^(KEbN0/10)/2)/sqrt(2);
%
    for i=1:Nsybl+Npath-1;
      Psig += abs(Y(i))^2/2;
      Pwgn += abs(Noise(i))^2;
    end
%
    Y = Y + Noise; %Add AWGN
    Yn = HH*Y; %Matched Filter
    Y2 = Yn(:,1);
    RdatI = zeros(1,Nsybl);
    RdatQ = zeros(1,Nsybl);
%Demodulation
    for i=1:Nsybl
      A = V(:,i);
      R = conj(dot(Y2,A));
      if real(R) > 0
        RdatI(i) = 1;
      elseif real(R) < 0
        RdatI(i) = 0;
      end
      if imag(R) > 0
        RdatQ(i) = 1;
      elseif imag(R) < 0
        RdatQ(i) = 0;
      end
    end
%Bit Error Rate
    for i=1:Nsybl
      if RdatI(i) == TdatI(i)
        Collect = Collect + 1;
      elseif RdatI(i) != TdatI(i)
        False = False + 1;
      end
      if RdatQ(i) == TdatQ(i)
        Collect = Collect + 1;
      elseif RdatQ(i) != TdatQ(i)
        False = False + 1;
      end
    end
  end
%
  EbN0 = 10*log10(Psig/Pwgn*2) %QPSK rate = 2 
  BER = False/(Collect + False)
  if BER > 0
    fprintf(fileID,'%f %f\n',EbN0,BER);  
  end
end
%
fclose(fileID);
toc;