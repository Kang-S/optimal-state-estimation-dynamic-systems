% This example uses an extended Kalman filter to estimate 
% the states of Van der Pol's. This program provides plots 
% of the simulated differenced measurements, the EKF 
% velocity estimates, and position and velocity errors 
% with 3-sigma outliers.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 3.5

% Written by John L. Crassidis 9/03

% Other Required Routines: polfun.m

% State and Initialize
dt=0.01;t=[0:dt:10]';m=length(t);
h=[1 0];r=0.01^2;
xe=zeros(m,2);x=zeros(m,2);p_cov=zeros(m,2);ym=zeros(m,1);
x0=[1;0];x(1,:)=x0';xe(1,:)=x0';
p0=1000*eye(2);p=p0;p_cov(1,:)=diag(p0)';

% True and Assumed Parameters
c=1;k=1;
cm=1.5;km=1.2;

% Process Noise (note: there is coupling but is ignored) 
q=0.2*[0 0;0 1];

% Main Routine
for i=1:m-1;

% Truth and Measurements
f1=dt*polfun(x(i,:),c,k);
f2=dt*polfun(x(i,:)+0.5*f1',c,k);
f3=dt*polfun(x(i,:)+0.5*f2',c,k);
f4=dt*polfun(x(i,:)+f3',c,k);
x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
ym(i)=x(i,1)+sqrt(r)*randn(1);

% Kalman Update
gain=p*h'*inv(h*p*h'+r);
p=(eye(2)-gain*h)*p;
xe(i,:)=xe(i,:)+gain'*(ym(i)-xe(i,1));

% Propagation
f1=dt*polfun(xe(i,:),cm,km);
f2=dt*polfun(xe(i,:)+0.5*f1',cm,km);
f3=dt*polfun(xe(i,:)+0.5*f2',cm,km);
f4=dt*polfun(xe(i,:)+f3',cm,km);
xe(i+1,:)=xe(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
fpart=[0 1;-4*cm*xe(i,1)*xe(i,2)-km -2*cm*(xe(i,1)^2-1)];
phi=c2d(fpart,[0;1],dt);
p=phi*p*phi'+q*dt;
p_cov(i+1,:)=diag(p)';

end

% 3-Sigma Outlier
sig3=p_cov.^(0.5)*3;

% Difference Measurement
ymd=diff(ym)/dt;
ymd(m)=ym(m-1);

% Plot Results
subplot(221)
plot(t,ymd)
set(gca,'Fontsize',12);
axis([0 10 -5 5]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-5 -2.5 0 2.5 5]);
xlabel('Time (Sec)')
ylabel('Differenced Measurement')

subplot(222)
plot(t,xe(:,2))
set(gca,'Fontsize',12);
axis([0 10 -5 5]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-5 -2.5 0 2.5 5]);
xlabel('Time (Sec)')
ylabel('Estimate')

subplot(223)
plot(t,xe(:,1)-x(:,1),t,sig3(:,1),t,-sig3(:,1))
set(gca,'Fontsize',12);
axis([0 10 -0.03 0.03]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-0.03 -0.02 -0.01 0 0.01 0.02 0.03]);
xlabel('Time (Sec)')
ylabel('Position Errors')

subplot(224)
plot(t,xe(:,2)-x(:,2),t,sig3(:,2),t,-sig3(:,2))
set(gca,'Fontsize',12);
axis([0 10 -0.5 0.5]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-0.5 -0.25 0 0.25 0.5]);
xlabel('Time (Sec)')
ylabel('Velocity Errors')