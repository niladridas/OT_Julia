clear all;close;clc
load a.dat; load b.dat; load K.dat; load lambda.dat; load U.dat;
% tic
t = cputime
[D,L,u,v] = sinkhornTransport(a,b,K,U,lambda,'marginalDifference',Inf,.5e-2,500,0);
% toc
e = cputime-t