function f = himmelblau_fun(alp,s,x);

% Written by John L. Crassidis 9/03

% Function Routine for Himmelblau

xx=x+alp*s;
f=(xx(1)^2+xx(2)-11)^2+(xx(1)+xx(2)^2-7)^2;