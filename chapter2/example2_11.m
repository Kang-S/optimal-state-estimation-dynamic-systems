% This example provides a plot of the ridge variance,
% the least squares variance, the ridge residual sum 
% squares, and the bias-squared quantities as a function 
% of the ridge parameter for a simple system. It outputs 
% the upper bound, the minimum value and the zero 
% derivative condition.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 2.11

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Ridge Variables
sig2=2;
lam=2;
x=1.5;
phi=[0:0.01:10]';
fac=(lam+phi);
p=lam*sig2./fac./fac;
q=(sig2*lam+phi.*phi*x*x)./fac./fac;
bias2=phi.*phi./fac./fac*x*x;
pls=sig2/lam;

% Plot Results
t=phi;
plot(t,q,t,ones(length(phi),1)*pls,t,p,'--',t,bias2,'-.')
axis([0 10 0 1.6]);
set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');
set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6])

text(5.5,1.05,'Least Squares Variance','Fontsize',12)
text(0.57,0.76,'Ridge','Fontsize',12)
text(4,0.15,'Ridge Variance','Fontsize',12)
text(2.3,0.6,'Bias-Squared','Fontsize',12)
h=get(gca,'Xlabel');
set(h,'String','\fontsize{12} {Ridge Parameter \phi}')
xlabel('Ridge Parameter')

% Upper Bound
upper_bound=2*sig2/(x*x-sig2/lam)

% Check Minimum
mink=find(q==min(q));
minphi=phi(mink)
facmin3=(lam+minphi)^3;
should_be_near_zero=2*minphi*lam*x^2/facmin3-2*sig2*lam/facmin3