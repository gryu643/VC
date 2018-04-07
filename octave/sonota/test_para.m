clear
pkg load parallel;
tic;

A = [1 2;3 4];
B = inv(A);

C = parallel_run(1:4,para(B,A))

toc;