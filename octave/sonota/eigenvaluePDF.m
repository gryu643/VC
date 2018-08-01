close all;
clear all;
tic;
Nsybl = 32;
Npath = 8;
PPLloop = 500;
trial = 1000;
stride=0.01;

s=0.0;
en=50.0;
result_rank=round((en-s)/stride)+1;
result=zeros(result_rank,2);

%Channel gain parameter
for i=1:Npath
  Ampd(i) = sqrt(1/Npath); %Equal Gain
%  Ampd(i) = sqrt(1/(2^(i-1))); %Exp. atten.
end

for loop=1:trial
    for i=1:Npath
        Cpath(i) = complex(randn,randn)/sqrt(2)*Ampd(i);
    end
%
    for i=1:Nsybl
        for j=1:Npath
            H(i+j-1,i) = Cpath(j);
        end
    end
% 
    for i=1:Nsybl+Npath-1
        for j = 1:Npath
            HE(i+j-1,i) = Cpath(j);
        end
    end
%
%    for l=1:Nsybl
%      for k = 1:Nsybl
%        Xppl(l,k) = 1.0 + 0.0*i;
%      end
%    end

%    [V,D] = PPL (H, HE, Xppl, Nsybl, Npath, PPLloop);
    
    HH = ctranspose(H);
    HHH = HH*H;
    [V,D] = eig(HHH);

    for i=1:Nsybl

        out = D(i,i);
        for j=1:result_rank
            aaa = s+stride*(j-1);
            bbb = s+stride*j;
            if (out>=aaa) && (out < bbb)
                result(j,2) = result(j,2) + 1.0;
            endif
        end
    end
end


for i=1:result_rank;
    result(i,1) = s + stride * (i-1);
    result(i,2) = result(i,2) /  (trial*32);
end

csvwrite('oct_evPDFeig.csv', result);

toc;