% This example shows a total least squares solution to fit the
% cooefficients of a polynomial function with noise on the measurements and
% basis functions

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 2.13

% Written by John L. Crassidis 10/09

% Other Required Routines: correlated_noise.m

% Time
t=[0:0.01:10]';m=length(t);

% True Output and Basis Functions
x_true=[1;0.5;0.3];
y=x_true(1)+x_true(2)*sin(t)+x_true(3)*cos(t);
h=[ones(m,1) sin(t) cos(t)];

% Measurement Covariance and Cholesky Decomposition
r=[1e-4 1e-6 1e-5 1e-9
    1e-6 1e-2 1e-7 1e-6
    1e-5 1e-7 1e-3 1e-6
    1e-9 1e-6 1e-6 1e-4];
c=chol(r);ci=inv(c);

% Number of Monte Carlo Runs and Storage of Variables
m_monte=1000;
x_itls_storage=zeros(m_monte,3);
x_tls_storage=zeros(m_monte,3);

% Main Loop
for i=1:m_monte

% Generate Noise
noise=correlated_noise(r,m);

ym=y+noise(:,4);
hm=[h(:,1)+noise(:,1) h(:,2)+noise(:,2) h(:,3)+noise(:,3)];

% Cholesky Decomposition and SVD of Augmented Matrix
d=[hm ym];
[m,n_d]=size(d);n=n_d-1;
d_star=d*ci;
[u,s,v]=svd(d_star,0);

% Isotropic Solution
x_itls=-v(1:n,n+1:n_d)*inv(v(n+1:n_d,n+1:n_d));

% Final Solution
x_tls=(ci(1:n,1:n)*x_itls-ci(1:n,n+1:n_d))*inv(ci(n+1:n_d,n+1:n_d));

% Store Solutions
x_itls_storage(i,:)=x_itls';
x_tls_storage(i,:)=x_tls';

end

% Show Results
mean_err_itls=mean(x_itls_storage)-x_true'
mean_err_tls=mean(x_tls_storage)-x_true'

cov_itls=cov(x_itls_storage);
cov_tls=cov(x_tls_storage);

trace_mse_itls=trace(cov_itls+mean_err_itls'*mean_err_itls)'
trace_mse_tls=trace(cov_tls+mean_err_tls'*mean_err_tls)'

% Compute Covariance
rot_cov=ci'*r*ci;
s_bar_inv=inv(s(1:n,1:n));
v_last_col_tran=v(:,n_d)';

d_cov=zeros(n_d,n_d);
for i = 1:m
 omega_mat=kron(v_last_col_tran,[s_bar_inv*u(i,1:n)';0]);
 d_cov=d_cov+omega_mat*rot_cov*omega_mat';
end
b_cov=v*d_cov*v';

cov_anal_itls=v(n_d,n_d)^(-2)*(b_cov(1:n,1:n)+v(n_d,n_d)^(-2)*v(1:n,n_d)*v(1:n,n_d)'*b_cov(n_d,n_d)-v(n_d,n_d)^(-1)*(b_cov(1:n,n_d)*v(1:n,n_d)'+v(1:n,n_d)*b_cov(1:n,n_d)'));
cov_anal_tls=ci(1:n,1:n)*cov_anal_itls*ci(1:n,1:n)'/ci(n_d,n_d)^2;

% Three Sigma Bounds
sig3=3*diag(cov_anal_tls).^(0.5);

% Plot Results
subplot(221)
plot(t,ym,t,hm(:,1),'--',t,hm(:,2),'-.',t,hm(:,3),':')
set(gca,'fontsize',12)
legend('y','h1','h2','h3',4)
ylabel('Measurements')
xlabel('Time (Sec)')

subplot(222)
x_axis=[1:m_monte]';
plot(x_axis,sig3(1)*ones(m_monte,1),'r--',x_axis,x_tls_storage(:,1)-x_true(1),'b',x_axis,-sig3(1)*ones(m_monte,1),'r--')
set(gca,'fontsize',12)
axis([0 m_monte -8e-3 8e-3])
ylabel('x1')
xlabel('Trial Run Number')

subplot(223)
plot(x_axis,sig3(2)*ones(m_monte,1),'r--',x_axis,x_tls_storage(:,2)-x_true(2),'b',x_axis,-sig3(2)*ones(m_monte,1),'r--')
set(gca,'fontsize',12)
axis([0 m_monte -0.01 0.01])
ylabel('x2')
xlabel('Trial Run Number')

subplot(224)
plot(x_axis,sig3(3)*ones(m_monte,1),'r--',x_axis,x_tls_storage(:,3)-x_true(3),'b',x_axis,-sig3(3)*ones(m_monte,1),'r--')
set(gca,'fontsize',12)
axis([0 m_monte -0.01 0.01])
ylabel('x3')
xlabel('Trial Run Number')