% This example uses an ensemble Kalman filter to estimate the
% states of a discretized one-dimensional diffusion model. 

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.7

% Written by John L. Crassidis 2/10

% Other Required Routines: none

% Define Parameters for Model and Time
length_y=4;delta_y=0.005;n=length_y/delta_y+1;
dt=0.01;tf=1;t=[0:dt:tf]';m=length(t);

% Get State Matrix
f=zeros(n);
f(1,1:2)=[-1 1];
f(n,n-1:n)=[1 -1];
for i=2:n-1,
 f(i,i-1:i+1)=[1 -2 1];
end
f=f/(delta_y)^2;

% Discrete-Time State Matrix
phi=c2d(f,zeros(n,1),dt);

% Process Noise Covariance
g=zeros(n,2);g(1,1)=1;g(n,2)=1;
q=1*eye(2);qd=dt*q;

% Initial Conditions
x0=1+[1:n]'*length_y/n;
x=zeros(m,n);x(1,:)=x0(:)';

% Output Matrix 
h=zeros(2,n);h(1,1)=1;h(2,n)=1;

% Measurements 
r=0.01*eye(2);
y=zeros(m,2);y(1,:)=(h*x(1,:)')';
ym=zeros(m,2);
ym(1,:)=y(1,:)+[sqrt(r(1,1))*randn(1) sqrt(r(2,2))*randn(1)];

% Number of Ensembles
n_ens=50;

% Initial Covariance and Ensemble Generation
p0=0.1^2*eye(n);
x_samp=kron(diag(p0).^(0.5),ones(1,n_ens)).*randn(n,n_ens)+kron(x0(:),ones(1,n_ens));
x_ens=x_samp;

% Estimates
xe=zeros(m,n);xe(1,:)=mean(x_ens,2)';
p_cov=zeros(m,n);

% Main Loop
for i =1:m-1
% Generate Truth
 x(i+1,:)=(phi*x(i,:)'+dt*g*[sqrt(qd(1,1))*randn(1);sqrt(qd(2,2))*randn(1)])';
 y(i+1,:)=(h*x(i+1,:)')';
 ym(i+1,:)=y(i+1,:)+ [sqrt(r(1,1))*randn(1) sqrt(r(2,2))*randn(1)];
 
% Ensemble Process Noise and Measurement Noise
 w_samp=kron(diag(qd).^(0.5),ones(1,n_ens)).*randn(length(qd),n_ens);
 v_samp=kron(diag(r).^(0.5),ones(1,n_ens)).*randn(length(r),n_ens);

% Compute Sample Covariances
 e_state=x_ens-kron(xe(i,:)',ones(1,n_ens));
 e_out=h*e_state;
 p_xy=1/(n_ens-1)*e_state*e_out';
 p_yy=1/(n_ens-1)*e_out*e_out';
 p_cov(i,:)=1/(n_ens-1)*diag(e_state*e_state')';

% Ensemble Kalman Update 
 gain_ens=p_xy*inv(p_yy);
 ens_res=kron(ym(i,:)',ones(1,n_ens))+v_samp-h*x_ens;
 x_ens=x_ens+gain_ens*ens_res;
 
% Ensemble Propagation and Estimate
 x_ens=phi*x_ens+g*w_samp;
 xe(i+1,:)=mean(x_ens,2)';
 
end

% Covariance at Final Point
e_state=x_ens-kron(xe(i+1,:)',ones(1,n_ens));
p_cov(i+1,:)=1/(n_ens-1)*diag(e_state*e_state')';
sig3=p_cov.^(0.5)*3;

% Plot Results
k_skip=[1:10:801]';
plot(t,xe(:,k_skip))
set(gca,'fontsize',12)
axis([0 1 1 5])
xlabel('Time (Sec)')
ylabel('Temperature')