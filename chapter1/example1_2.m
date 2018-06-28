% This example illustrates the basic concept of using linear 
% least squares to estimate the parameters of a simple dynamic
% system. The program provides a plot of the measurements used 
% in least squares and the best fit.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 1.2

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Get Truth and Measurements
dt=0.1;tf=10;
t=[0:dt:tf];
m=length(t);
y=zeros(m,1);
a=-1;b=1;c=1;d=0;u=[100;zeros(m-1,1)];
[ad,bd]=c2d(a,b,dt);
y=dlsim(ad,bd,c,d,u);
ym=y+0.08*randn(m,1);

% Least Squares Estimate
w=inv(0.08^2);
h=[ym(1:m-1) u(1:m-1)];
xe=inv(h'*w*h)*h'*w*ym(2:m)

% Plot Results
ye=dlsim(xe(1),xe(2),c,d,u);
plot(t,ym,'*',t,ye)
set(gca,'fontsize',12);
xlabel('Measurements of {\it y}({\it t}) and Best Fit')
legend('Measurements','Best Fit')
set(gca,'Ytick',[-2 0 2 4 6 8 10]);
set(gca,'Xtick',[0 1 2 3 4 5 6 7 8 9 10]);