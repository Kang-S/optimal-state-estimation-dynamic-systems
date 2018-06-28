% This example uses a SVD decomposition to solve a 
% simple inequality constrained problem. It provides 
% the unconstrained least squares solution for comparison 
% purposes.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 1.9

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Measurements and Other Parameters
t=[0:0.1:10]';
y=3+2*t+t.*t;
ym=y+randn(length(t),1);
h=[ones(length(t),1) t t.*t];
[u,s,v]=svd(h,0);
z=u'*ym;
r=rank(h);

le=0;
for i=1:length(z)
 le=le+(z(i)/s(i,i))^2;
end

% Value for gamma
gam=sqrt(1^2+2^2+3^2);

% Main SVD Algorithm
xe=zeros(3,1);
lam=0;fc=0;fdc=0;fstop=100;lami=0;
if le > gam*gam
 while(norm(fstop)>1e-8)
  lami=lami+1;
  for i=1:length(z);
   fc=fc+(s(i,i)*z(i)/(s(i,i)^2+lam))^2;
   fdc=fdc-2*(s(i,i)*z(i))^2*(s(i,i)^2+lam)^(-3);
  end
 
  fc=fc-gam*gam;
  fstop=fc/fdc;
  lam=lam-fc/fdc;

  if (lami>20000) break; end
 end

 for i=1:length(z);
  xe=xe+s(i,i)*z(i)/(s(i,i)^2+lam)*v(:,i);
 end

else

 for i=1:length(z)
  xe=xe+(z(i)/s(i,i))*v(:,i);
 end

end

xe_svd=xe
prod_svd=xe_svd'*xe_svd
disp(' ')

% Least Square Solution
xe_ls=inv(h'*h)*h'*ym
prod_ls=xe_ls'*xe_ls