function f=athansfun(x,gam)

% Written by John L. Crassidis 9/03

% Function Routine for Athans Problem
[m,n]=size(x);
f=zeros(m,n);
f(1,:)=-x(2,:);
f(2,:)=-exp(-gam*x(1,:)).*(x(2,:).^2).*x(3,:);
f(3,:)=zeros(1,n);