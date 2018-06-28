function [rho,limit95]=whiteness_test(y,point)

% This program gives statistics for a whiteness test.
% The inputs are:
%        y = signal 
%    point = desired point to test
% The outputs are:
%      rho = autocorrelation
%  limit95 = 95% limit 

% Written by John L. Crassidis 9/03

[m,n]=size(y);

bigk=m-point+1;

dot1=sum(dotprod1(y(1:bigk,:),y(1:bigk,:)));
dot2=sum(dotprod1(y(point:m,:),y(point:m,:)));
dot3=sum(dotprod1(y(1:bigk,:),y(point:m,:)));

rho=abs(1/sqrt(n)*dot3*(dot1)^(-0.5)*(dot2)^(-0.5));

limit95=1.96*1/sqrt(bigk);