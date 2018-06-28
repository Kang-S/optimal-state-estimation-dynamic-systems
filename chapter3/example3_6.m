% This example uses an extended Kalman filter to estimate 
% the damping coefficient in Van der Pol's. This program 
% provides plots of the velocity estimates for various values 
% of the process noise variance, and parameter estimate errors
% with 3-sigma outliers.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 3.6

% Written by John L. Crassidis 9/03

% Other Required Routines: polfun.m, polfun_id.m

% State and Initialize
dt=0.01;t=[0:dt:10]';m=length(t);
h=[1 0 0];r=0.01^2;
x=zeros(m,2);ym=zeros(m,1);
h=[1 0 0];
xe_id=zeros(m,3);p_cov_id=zeros(m,3);
x0=[1;0];x(1,:)=x0';xe1(1,:)=[1;0;0]';xe2(1,:)=[1;0;0]';
p0=1000*eye(3);p=p0;p_cov1(1,:)=diag(p0)';p_cov2(1,:)=diag(p0)';

% True and Assumed Parameters
c=1;k=1;km=k;

% Truth and Measurements
for i=1:m-1;

f1=dt*polfun(x(i,:),c,k);
f2=dt*polfun(x(i,:)+0.5*f1',c,k);
f3=dt*polfun(x(i,:)+0.5*f2',c,k);
f4=dt*polfun(x(i,:)+f3',c,k);
x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
ym(i)=x(i,1)+sqrt(r)*randn(1);

end

% Run for Different Values of q
for j = 1:2

if j == 1,
% Main Routine
for i=1:m-1;

% Kalman Update
gain=p*h'*inv(h*p*h'+r);
p=(eye(3)-gain*h)*p;
xe1(i,:)=xe1(i,:)+gain'*(ym(i)-xe1(i,1));

% Propagation
f1=dt*polfun_id(xe1(i,:),km);
f2=dt*polfun_id(xe1(i,:)+0.5*f1',km);
f3=dt*polfun_id(xe1(i,:)+0.5*f2',km);
f4=dt*polfun_id(xe1(i,:)+f3',km);
xe1(i+1,:)=xe1(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
fpart=[0 1 0
-4*xe1(i,3)*xe1(i,1)*xe1(i,2)-k -2*xe1(i,3)*(xe1(i,1)^2-1) -2*(xe1(i,1)^2-1)*xe1(i,2)
0 0 0];
phi=c2d(fpart,[0;0;1],dt);
q=.1*[0 0 0;0 0 0;0 0 1];
p=phi*p*phi'+q*dt;
p_cov1(i+1,:)=diag(p)';

end

% 3-Sigma Outlier
sig3_id1=(p_cov1.^(0.5))*3;

else

p=p0;    
% Main Routine
for i=1:m-1;

% Kalman Update
gain=p*h'*inv(h*p*h'+r);
p=(eye(3)-gain*h)*p;
xe2(i,:)=xe2(i,:)+gain'*(ym(i)-xe2(i,1));

% Propagation
f1=dt*polfun_id(xe2(i,:),km);
f2=dt*polfun_id(xe2(i,:)+0.5*f1',km);
f3=dt*polfun_id(xe2(i,:)+0.5*f2',km);
f4=dt*polfun_id(xe2(i,:)+f3',km);
xe2(i+1,:)=xe2(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
fpart=[0 1 0
-4*xe2(i,3)*xe2(i,1)*xe2(i,2)-k -2*xe2(i,3)*(xe2(i,1)^2-1) -2*(xe2(i,1)^2-1)*xe2(i,2)
0 0 0];
phi=c2d(fpart,[0;0;1],dt);
q=.01*[0 0 0;0 0 0;0 0 1];
p=phi*p*phi'+q*dt;
p_cov2(i+1,:)=diag(p)';

end

sig3_id2=(p_cov2.^(0.5))*3;

end

end

% Plot Results
subplot(221)
plot(t,xe1(:,2),'--',t,x(:,2))
set(gca,'Fontsize',12);
axis([0 10 -5 5]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-5 -2.5 0 2.5 5]);
xlabel('Time (Sec)')
ylabel('Velocity Estimate')
text(4,-3,'{\it q} = 0.1')

subplot(222)
plot(t,xe2(:,2),'--',t,x(:,2))
set(gca,'Fontsize',12);
axis([0 10 -5 5]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-5 -2.5 0 2.5 5]);
xlabel('Time (Sec)')
ylabel('Velocity Estimate')
text(4,-3,'{\it q} = 0.01')

subplot(223)
plot(t,xe1(:,3)-1,t,sig3_id1(:,3),t,-sig3_id1(:,3))
set(gca,'Fontsize',12);
axis([0 10 -3 3]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-3 -2 -1 0 1 2 3]);
xlabel('Time (Sec)')
ylabel('Parameter Error')
text(4,-2,'{\it q} = 0.1')

subplot(224)
plot(t,xe2(:,3)-1,t,sig3_id2(:,3),t,-sig3_id2(:,3))
set(gca,'Fontsize',12);
axis([0 10 -3 3]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-3 -2 -1 0 1 2 3]);
xlabel('Time (Sec)')
ylabel('Parameter Error')
text(4,-2,'{\it q} = 0.01')