function fun=ric_backfun(t,x);

% Written by John L. Crassidis 9/03

% Function Routine for Backward Riccati Equation

f=-1;q=2;r=1;
fun=2*f*x-x^2*q+inv(r);