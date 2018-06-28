% This example illustrates the basic concept of using linear 
% least squares for curve fitting a set of measurements. The program 
% provides a plot of the measurements used in least squares and the best fit.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins 
% Example 1.1

% Written by John L. Crassidis 2/10

% Other Required Routines: none

% Generate Measurements
t=[0:0.1:10]';m=length(t);
r=0.001;
ym=0.3*sin(t)+0.5*cos(t)+0.1*t+sqrt(r)*randn(m,1);
 
% Form Basis Function Matrix
h=[sin(t) cos(t) t cos(t).*sin(t) t.^2];
 
% Least Squares Estimate
xe=inv(h'*h)*h'*ym
 
% Estimated Output
ye=xe(1)*sin(t)+xe(2)*cos(t)+xe(3)*t+xe(4)*cos(t).*sin(t)+xe(5)*t.^2;
 
% Plot Results
plot(t,ym,'*',t,ye)
set(gca,'fontsize',12)
legend('Measurements','Best Fit','Location','SouthEast')
xlabel('Time (Sec)')
ylabel('Measurements and Best Fit')