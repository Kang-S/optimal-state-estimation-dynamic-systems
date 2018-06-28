% This example uses the parameters from example 1.2, but 
% solves the problem using sequential least squares. The 
% program provides plots of the estimates and covariances.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 1.5

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% Get Truth and Measurements
dt=0.1;tf=10;
t=[0:dt:tf];
m=length(t);
a=-1;b=1;c=1;d=0;u=[100;zeros(m-1,1)];
[ad,bd]=c2d(a,b,dt);
y=dlsim(ad,bd,c,d,u);
ym=y+0.08*randn(m,1);

% Weight and H Matrix
w=inv(0.08^2);
h=[ym(1:m-1) u(1:m-1)];

% Initial Conditions for Sequential Algorithm
alpha=1e3;
beta=[1e-2;1e-2];
p0=inv(1/alpha/alpha*eye(2)+h(1,:)'*w*h(1,:));
x0=p0*(1/alpha*beta+h(1,:)'*w*ym(2));

% Sequential Least Squares
xr=zeros(m-1,2);xr(1,:)=x0';
p=zeros(m-1,2);p(1,:)=diag(p0)';pp=p0;
for i=1:m-2;
 k=pp*h(i+1,:)'*inv(h(i+1,:)*pp*h(i+1,:)'+inv(w));
 pp=(eye(2)-k*h(i+1,:))*pp;
 xr(i+1,:)=xr(i,:)+(k*(ym(i+2)-h(i+1,:)*xr(i,:)'))';
 p(i+1,:)=diag(pp)';
end

% Plot Results
subplot(221)
plot(t(1:100),xr(:,1));grid
axis([0 10 0 10]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 2 4 6 8 10]);
set(gca,'ytick',[0 2 4 6 8 10]);
xlabel('Time (Sec)')
ylabel('{\it x}_1 Estimate')

subplot(222)
plot(t(1:100),xr(:,2));grid
axis([0 10 0.0945 0.096]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 2 4 6 8 10]);
set(gca,'ytick',[0.0945 0.095 0.0955 0.096]);
xlabel('Time (Sec)')
ylabel('{\it x}_2 Estimate')

subplot(223)
semilogy(t(1:100),p(:,1));grid
axis([0 10 1e-5 1e5])
set(gca,'xtick',[0 2 4 6 8 10]);
set(gca,'ytick',[1e-5 1 1e5])
xlabel('Time (Sec)')
ylabel('{\it P}_{11}')

subplot(224)
semilogy(t(1:100),p(:,2));grid
set(gca,'xtick',[0 2 4 6 8 10]);
xlabel('Time (Sec)')
ylabel('{\it P}_{22}')