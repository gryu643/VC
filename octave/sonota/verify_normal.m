close all;
clear all;
tic;

trial=1000000;
stride=0.01;
avg=0.0;
variance=0.0;
s=-5.0;
en=5.0;
result_rank=round((en-s)/stride)+1;
result=zeros(result_rank,2);

fileID = fopen('octave_normal.csv','w');

%implementation
out = rand;
for i=1:trial
    out = randn;
    for j=1:result_rank
        aaa = s+stride*(j-1);
        bbb = s+stride*j;
        if (out>=aaa) && (out < bbb)
            result(j,2) = result(j,2) + 1.0;
        endif
    end
end

for i=1:result_rank;
    result(i,1) = s + stride * (i-1);
    result(i,2) = result(i,2) / trial;
end

i=1:result_rank;
%fprintf(fileID,'%f, %f\n',result(i,1),result(i,2));
csvwrite('octave_normal.csv', result);

%calculate avg
for i=1:result_rank
    avg += result(i,1)*result(i,2);
end

%calculate variance
for j=1:result_rank
    variance += ((result(j,1)-avg)^2).*result(j,2);
end

avg
variance

%file close
fclose(fileID);
toc;