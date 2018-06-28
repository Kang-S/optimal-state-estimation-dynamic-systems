function f=f8_fun(x,c)

% Written by John L. Crassidis 9/03

% Main Function Routine for f8 Aircraft

% Coefficients
c1=c(1);c2=c(2);c3=c(3);c4=c(4);c5=c(5);c6=c(6);c7=c(7);
c8=c(8);c9=c(9);c10=c(10);c11=c(11);

% Function 
f=zeros(3,1);
f(1)=c1*x(3)-c2*x(1)^2*x(3)-c3*x(1)*x(3)-c4*x(1)+c5*x(1)^2+c6*x(1)^3-c7*x(2)^2;
f(2)=x(3);
f(3)=-c8*x(3)-c9*x(1)-c10*x(1)^2-c11*x(1)^3;