function xout=air_covfuns(x,f,q);

% Written by John L. Crassidis 9/03

% Function Routine for Covariance

x=x(:);
p=[x(1:7) x(8:14) x(15:21) x(22:28) x(29:35) x(36:42) x(43:49)];
pdot=f*p+p*f'+q;
xout=[pdot(:,1);pdot(:,2);pdot(:,3);pdot(:,4);pdot(:,5);pdot(:,6);pdot(:,7)];