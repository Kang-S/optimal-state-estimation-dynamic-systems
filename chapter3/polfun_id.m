function f=polfun_id(x,k)

% Written by John L. Crassidis 9/03

% Function Routine for Van der Pol's Equation (id purposes)

f=zeros(3,1);
f(1)=x(2);
f(2)=-2*x(3)*(x(1)^2-1)*x(2)-k*x(1);
f(3)=0;