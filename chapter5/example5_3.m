% This example uses a nonlinear RTS smoother to estimate 
% the states of Van der Pol's. This program provides plots 
% of the EKF and smoother estimates, and the EKF and smoother
% velocity errors with 3-sigma outliers.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 5.3

% Written by John L. Crassidis 9/03

% Other Required Routines: polfun.m, polfun_cov.m, polfun_b, polfun_covb

% State and Initialize
dt=0.01;t=[0:dt:10]';m=length(t);
h=[1 0];r=0.01^2;
xe=zeros(m,2);x=zeros(m,2);p_cov=zeros(m,2);p_cov_s=zeros(m,2);
ym=zeros(m,1);
x0=[1;0];x(1,:)=x0';xe(1,:)=x0';
p0=1000*eye(2);p=p0;p_cov(1,:)=diag(p0)';
p_storep=zeros(m,4); p_storeu=zeros(m,4);
xs=zeros(m,2);

% Truth and Model Parameters
c=1;k=1;
cm=1.5;km=1.2;

% Process Noise (note: now continuous)
q=.2*[0 0;0 1];

% Main Forward Routine
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
p_storeu(i,:)=[p(1,1) p(1,2) p(2,1) p(2,2)];
xe(i,:)=xe(i,:)+gain'*(ym(i)-xe(i,1));

% State Propagation
f1=dt*polfun(xe(i,:),cm,km);
f2=dt*polfun(xe(i,:)+0.5*f1',cm,km);
f3=dt*polfun(xe(i,:)+0.5*f2',cm,km);
f4=dt*polfun(xe(i,:)+f3',cm,km);
xe(i+1,:)=xe(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');

% Covariance Propagation
xdum=[p(1,1) p(1,2) p(2,1) p(2,2)];
f1=dt*polfun_cov(xdum,xe(i,:),cm,km,q);
f2=dt*polfun_cov(xdum+0.5*f1',xe(i,:),cm,km,q);
f3=dt*polfun_cov(xdum+0.5*f2',xe(i,:),cm,km,q);
f4=dt*polfun_cov(xdum+f3',xe(i,:),cm,km,q);
xdum=xdum+1/6*(f1'+2*f2'+2*f3'+f4');
p=[xdum(1) xdum(2);xdum(3) xdum(4)];

p_storep(i+1,:)=[p(1,1) p(1,2) p(2,1) p(2,2)];
p_cov(i+1,:)=diag(p)';

end

% RTS Initialize 
xs(m,:)=xe(m,:);
p_cov_s(m,:)=p_cov(m,:);
pb=[p_storep(m,1) p_storep(m,2) 
    p_storep(m,3) p_storep(m,4)];

p_prop=[p_storep(m,1) p_storep(m,2) 
        p_storep(m,3) p_storep(m,4)];
ddd_i=inv(h*p_prop*h'+r);
gain=p_prop*h'*ddd_i;
lam=-(h'*ddd_i*(ym(m)-xe(m,1)))';

% Main Backward Routine
for i=m-1:-1:1

% Covariance
p_prop=[p_storep(i+1,1) p_storep(i+1,2) 
p_storep(i+1,3) p_storep(i+1,4)];
p_propi=inv(p_prop);

% Backward State
f1=dt*polfun_b(xs(i+1,:),p_prop,xe(i+1,:),cm,km,q);
f2=dt*polfun_b(xs(i+1,:)+0.5*f1',p_prop,xe(i+1,:),cm,km,q);
f3=dt*polfun_b(xs(i+1,:)+0.5*f2',p_prop,xe(i+1,:),cm,km,q);
f4=dt*polfun_b(xs(i+1,:)+f3',p_prop,xe(i+1,:),cm,km,q);
xs(i,:)=xs(i+1,:)+1/6*(f1'+2*f2'+2*f3'+f4');

% Backward Covariance
xdum=[pb(1,1) pb(1,2) pb(2,1) pb(2,2)];
f1=dt*polfun_covb(xdum,xe(i+1,:),cm,km,q,p_propi);
f2=dt*polfun_covb(xdum+0.5*f1',xe(i+1,:),cm,km,q,p_propi);
f3=dt*polfun_covb(xdum+0.5*f2',xe(i+1,:),cm,km,q,p_propi);
f4=dt*polfun_covb(xdum+f3',xe(i+1,:),cm,km,q,p_propi);
xdum=xdum+1/6*(f1'+2*f2'+2*f3'+f4');
pb=[xdum(1) xdum(2);xdum(3) xdum(4)];
p_cov_s(i,:)=diag(pb)';

end

% 3-Sigma Outliers
sig3=p_cov.^(0.5)*3;
sig3_s=p_cov_s.^(0.5)*3;

% Plot Results
subplot(221)
plot(t,xe(:,2))
set(gca,'Fontsize',12);
axis([0 10 -5 5]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-5 -2.5 0 2.5 5]);
xlabel('Time (Sec)')
ylabel('EKF Estimate')

subplot(222)
plot(t,xs(:,2))
set(gca,'Fontsize',12);
axis([0 10 -5 5]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-5 -2.5 0 2.5 5]);
xlabel('Time (Sec)')
ylabel('Smoother Estimate')

subplot(223)
plot(t,xe(:,2)-x(:,2),t,sig3(:,2),t,-sig3(:,2))
set(gca,'Fontsize',12);
axis([0 10 -0.5 0.5]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-0.5 -0.25 0 0.25 0.5]);
xlabel('Time (Sec)')
ylabel('EKF Velocity Error')

subplot(224)
plot(t,xs(:,2)-x(:,2),t,sig3_s(:,2),t,-sig3_s(:,2))
set(gca,'Fontsize',12);
axis([0 10 -0.5 0.5]);
set(gca,'Xtick',[0 2 4 6 8 10]);
set(gca,'Ytick',[-0.5 -0.25 0 0.25 0.5]);
xlabel('Time (Sec)')
ylabel('Smoother Velocity Error')