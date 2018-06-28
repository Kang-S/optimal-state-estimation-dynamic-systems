function f = quadratic_fun(alp,s,x);

% Written by John L. Crassidis 9/03

% Function Routine for Quadratic

xx=x+alp*s;
f=4*xx(1)^2+3*xx(2)^2-4*xx(1)*xx(2)+xx(1);