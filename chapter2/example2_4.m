% This example solves a simple constrained least squares problem. 
% This program provides plots of the parameter estimate errors and 
% 3-sigma outliers.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 2.4

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Truth
t=[0:6/90:2]';m=length(t);
y=t+sin(t)+2*cos(2*t);

% Pre-allocate Space
xe=zeros(1000,3);
pcov=zeros(1000,3);

% Monte Carlo Simulation
for i=1:1000,

% Measurements with Noise 
 sig=0.1;
 ym=y+sig*randn(m,1);   
 ym(1:3)=y(1:3);
 ym1=ym(4:31);
 ym2=ym(1:3);

% Basis Functions
 h1=[t(4:31) sin(t(4:31)) cos(2*t(4:31))];
 h2=[t(1:3) sin(t(1:3)) cos(2*t(1:3))];

% Get Constrained Estimate
 p_bar=inv(h1'*h1);
 xe_bar=p_bar*h1'*ym1;
 gain=p_bar*h2'*inv(h2*p_bar*h2');
 xe(i,:)=(xe_bar+gain*(ym2-h2*xe_bar))';   
 p=sig^2*(eye(3)-gain*h2)*p_bar*(eye(3)-gain*h2)';
 pcov(i,:)=diag(p)';
end

% Plot Results
subplot(221)
plot(t,ym,'*');
axis([0 2 0 3]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 0.5 1 1.5 2]);
set(gca,'ytick',[0 1 2 3]);
xlabel('Time (Sec)')
ylabel('Measurements')

subplot(222)
plot([1:1000],xe(:,1)-1,[1:1000],pcov(:,1).^(0.5)*3,[1:1000],-pcov(:,1).^(0.5)*3);
axis([0 1000 -1e-9 1e-9]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 250 500 750 1000]);
set(gca,'ytick',[-1e-9 -0.5e-9 0 0.5e-9 1e-9]);
xlabel('Trial Run Number')
ylabel('{\it x}_1 Error and 3{\sigma} Outlier')

subplot(223)
plot([1:1000],xe(:,2)-1,[1:1000],pcov(:,2).^(0.5)*3,[1:1000],-pcov(:,2).^(0.5)*3);
axis([0 1000 -1e-9 1e-9]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 250 500 750 1000]);
set(gca,'ytick',[-1e-9 -0.5e-9 0 0.5e-9 1e-9]);
xlabel('Trial Run Number')
ylabel('{\it x}_2 Error and 3{\sigma} Outlier')

subplot(224)
plot([1:1000],xe(:,3)-2,[1:1000],pcov(:,3).^(0.5)*3,[1:1000],-pcov(:,3).^(0.5)*3);
axis([0 1000 -1e-9 1e-9]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 250 500 750 1000]);
set(gca,'ytick',[-1e-9 -0.5e-9 0 0.5e-9 1e-9]);
xlabel('Trial Run Number')
ylabel('{\it x}_3 Error and 3{\sigma} Outlier')