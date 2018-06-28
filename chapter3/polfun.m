function f=polfun(x,c,k)

% Written by John L. Crassidis 9/03

% Function Routine for Van der Pol's Equation

f=[x(2);-2*c*(x(1)^2-1)*x(2)-k*x(1)];