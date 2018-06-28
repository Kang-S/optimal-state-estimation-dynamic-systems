% This program runs an MMAE estimator to determine the 
% process noise variance of a linear system. It outputs the 
% parameter estimate and MMAE overal state estimate

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.5

% Written by John L. Crassidis 1/10

% Other Required Routines: none

% Initialize
dt=0.01;tf=10;t=[0:dt:tf]';m=length(t);
x0=[1;1];

% System Matrices
a=[0 1;-3 -3];b=[0;1];h=[1 0];
r=0.01;q=10;
phi=c2d(a,b,dt);
upsilon=[0;0.01];
n=length(phi);

% Generate Measurements
x=zeros(m,2);x(1,:)=x0';
y=zeros(m,1);y(1)=h*x(1,:)';
for i = 1:m-1
 x(i+1,:)=(phi*x(i,:)'+upsilon*sqrt(q)*randn(1))';
 y(i+1)=h*x(i,:)';
end
ym=y+sqrt(r)*randn(m,1);

% Trial Values
q_trial=[0.1;0.5;1;10;20;100;1000;1e4;1e5]';
m_trial=length(q_trial);

% Weights and Initial Estimates for Parameters
weights=ones(m_trial,1)/m_trial;
param_est=zeros(m,1);

% Initial Covariance 
p=diag([0.001^2 0.001^2]);

% MMAE State Initialization
xe=zeros(m,n);xe(1,:)=x0';

% MMAE Storage of State and Covariance
xx_prop=kron(xe(1,:)',ones(1,m_trial));xx_up=zeros(n,m_trial);
pp_prop=kron(ones(1,m_trial),p);pp_up=zeros(n,n*m_trial);
like=zeros(m_trial,1);

% Allocate MMAE State and Parameter Covariances
pcov_x=zeros(m,n);pcov_p=zeros(m,1);

% Main Loop 
for i=1:m
    
 for j_trial=1:m_trial

% Likelihood
  res=ym(i)-h*xx_prop(:,j_trial);
  e_cov=h*pp_prop(:,n*j_trial-n+1:n*j_trial)*h'+r;
  like(j_trial)=1/sqrt(det(2*pi*e_cov))*exp(-0.5*res'*inv(e_cov)*res);    
     
% Kalman Gain    
  gain=pp_prop(:,n*j_trial-n+1:n*j_trial)*h'*inv(h*pp_prop(:,n*j_trial-n+1:n*j_trial)*h'+r);

% Update
  xx_up(:,j_trial)=xx_prop(:,j_trial)+gain*(ym(i)-h*xx_prop(:,j_trial));
  pp_up(:,n*j_trial-n+1:n*j_trial)=[eye(n)-gain*h]*pp_prop(:,n*j_trial-n+1:n*j_trial);

% Propagate
  xx_prop(:,j_trial)=phi*xx_up(:,j_trial);
  pp_prop(:,n*j_trial-n+1:n*j_trial)=phi*pp_up(:,n*j_trial-n+1:n*j_trial)*phi'+upsilon*q_trial(j_trial)*upsilon';

 end
 
% Updated Weights
 weights=weights.*like;
 weights=weights/sum(weights);

% MMAE State and Parameter Estimates
 xe(i,:)=sum(xx_up.*kron(weights',ones(n,1)),2)';
 param_est(i,:)=sum(weights'.*q_trial);

% MMAE State and Parameter Covariances
 x_diff=xx_up-kron(xe(i,:)',ones(1,m_trial));
 px=kron(weights',ones(n,n)).*pp_up;
 pmat_x=(x_diff.*kron(weights',ones(n,1)))*x_diff'+[sum(px(:,1:n:end),2) sum(px(:,2:n:end),2)];
 pcov_x(i,:)=diag(pmat_x)';
 
 p_diff=q_trial-param_est(i);
 pmat_p=(p_diff.*kron(weights',ones(1,1)))*p_diff';
 pcov_p(i,:)=diag(pmat_p)';
 
end

% Three-Sigam Bounds
sig3_x=pcov_x.^(0.5)*3;
sig3_p=pcov_p.^(0.5)*3;

% Plot Results
plot(t,sig3_x(:,1),'r',t,xe(:,1)-x(:,1),'b',t,-sig3_x(:,1),'r')
set(gca,'fontsize',12)
xlabel('Time (Sec)')
ylabel('x1 Errors')

pause

plot(t,sig3_x(:,2),'r',t,xe(:,2)-x(:,2),'b',t,-sig3_x(:,2),'r')
set(gca,'fontsize',12)
xlabel('Time (Sec)')
ylabel('x2 Errors')

pause

plot(t,sig3_p,'r',t,param_est-q,'b',t,-sig3_p,'r')
axis([0 10 -100 100])
set(gca,'fontsize',12)
xlabel('Time (Sec)')
ylabel('Parameter Estimate Errors')