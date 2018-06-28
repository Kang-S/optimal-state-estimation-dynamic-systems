function f=orbitfun(x,mu);

% Written by John L. Crassidis 9/03

% Function Routine for Orbital Equations

f=zeros(6,1);
r=sqrt(x(1)^2+x(2)^2+x(3)^2);
r3=r^3;
f(1)=x(4);
f(2)=x(5);
f(3)=x(6);
f(4)=-mu/r3*x(1);
f(5)=-mu/r3*x(2);
f(6)=-mu/r3*x(3);