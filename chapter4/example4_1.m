% This example shows the design results of a simple colored-noise 
% filter using the longitudinal short-period dynamics of an aircraft. 
% This program provides a plot of the 3-sigma outliers for various 
% values of the process noise covariance and gust noise parameter.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.1

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Variance, Natural Frequency and Damping Ratio
r=(1*pi/180)^2;
omega=1;
zeta=sqrt(2)/2;

% Gust and Process Noise Parameters
q=[1:2:50]';
a=[0:1:10]';
out=zeros(length(a),length(q));

% Main Loop
for i=1:length(a)
 for j =1:length(q)  
  f=[0 1 0;-omega^2 -2*zeta*omega 1;0 0 -a(i)];
  g=[0;0;1];
  h=[1 0 0];
  p=are(f',h'*inv(r)*h,g*q(j)*g'); 
  out(i,j)=p(1,1)^(0.5)*3*180/pi;
 end
end

% Plot Results
mesh(a,q,out')
axis([0 10 0 50 4 11]);
set(gca,'Fontsize',12);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[0 10 20 30 40 50]);
set(gca,'Ztick',[4 5 6 7 8 9 10 11]);
xlabel('Gust Noise Parameter {\it a}')
ylabel('Covariance {\it q}')
zlabel('3{\sigma} Pitch Outlier (Deg)')