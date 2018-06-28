% This example solves a simple least squares problem. This program
% provides plots of the parameter estimate errors and 3-sigma 
% outliers.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 2.2

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% True System
dt=0.01;tf=10;
t=[0:dt:tf]';
m=length(t);
y=cos(t)+2*sin(t)+cos(2*t)+2*sin(3*t);

% Pre-allocate Space
xe=zeros(1000,4);
pcov=zeros(1000,4);

% Monte Carlo Simulation
for i=1:1000,
 ym=y+0.1*randn(m,1);w=1/0.01;
 h=[cos(t) sin(t) cos(2*t) sin(3*t)];
 p=inv(h'*w*h);
 xe(i,:)=(p*h'*w*ym)';
 pcov(i,:)=diag(p)';
end

% Plot Results
subplot(221)
plot([1:1000],xe(:,1)-1,[1:1000],pcov(:,1).^(0.5)*3,[1:1000],-pcov(:,1).^(0.5)*3);
axis([0 1000 -0.02 0.02]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 250 500 750 1000]);
set(gca,'ytick',[-0.02 -0.01 0 0.01 0.02]);
xlabel('Trial Run Number')
ylabel('{\it x}_1 Error and 3{\sigma} Outlier')

subplot(222)
plot([1:1000],xe(:,2)-2,[1:1000],pcov(:,2).^(0.5)*3,[1:1000],-pcov(:,2).^(0.5)*3);
axis([0 1000 -0.02 0.02]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 250 500 750 1000]);
set(gca,'ytick',[-0.02 -0.01 0 0.01 0.02]);
xlabel('Trial Run Number')
ylabel('{\it x}_2 Error and 3{\sigma} Outlier')

subplot(223)
plot([1:1000],xe(:,3)-1,[1:1000],pcov(:,3).^(0.5)*3,[1:1000],-pcov(:,3).^(0.5)*3);
axis([0 1000 -0.02 0.02]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 250 500 750 1000]);
set(gca,'ytick',[-0.02 -0.01 0 0.01 0.02]);
xlabel('Trial Run Number')
ylabel('{\it x}_3 Error and 3{\sigma} Outlier')

subplot(224)
plot([1:1000],xe(:,4)-2,[1:1000],pcov(:,4).^(0.5)*3,[1:1000],-pcov(:,4).^(0.5)*3);
axis([0 1000 -0.02 0.02]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 250 500 750 1000]);
set(gca,'ytick',[-0.02 -0.01 0 0.01 0.02]);
xlabel('Trial Run Number')
ylabel('{\it x}_4 Error and 3{\sigma} Outlier')
