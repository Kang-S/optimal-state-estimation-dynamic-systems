function f=polfun_cov(x,xe,c,k,q)

% Written by John L. Crassidis 9/03

% Function Routine for Van der Pol's Covariance Equation in Forward Integration

% State Matrix
fpart=[0 1;-4*c*xe(1)*xe(2)-k -2*c*(xe(1)^2-1)];

% Covariance
p=[x(1) x(2);x(3) x(4)];
pdot=fpart*p+p*fpart'+q;
f=[pdot(1,1);pdot(1,2);pdot(2,1);pdot(2,2)];