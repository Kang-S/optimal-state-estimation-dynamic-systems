% In this example the Rao-Blackwellized Particle Filter (RBPF) is used 
% to estimate the states of a finite impulse response (FIR) filter

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.13

% Written by John L. Crassidis 2/10

% Other Required Routines: resample_pf

% Time
dt=1;tf=100;t=[0:dt:tf]';m=length(t);

% Process Noise and Measurement Noise Variances
q1=0.09;
q2=0.04;
r=1e-3;

% Initial Condition
x=zeros(m,2);x0=[1;2];x(1,:)=x0(:)';
ym=zeros(m,1);ym(1)=x(1,1)*x(1,2)+sqrt(r)*randn(1);

% Particles
n=500;
xe=zeros(m,2);xe(1,:)=x0(:)';
w_particle=inv(n)*ones(n,1);

% Kalman Part
p0=1;p=p0*ones(n,1);
x_particle1=(cos(x0(1))+sin(x0(1)))+sqrt(q1)*randn(n,1);
x_particle2=sqrt(p0)*randn(n,1);
p_cov=zeros(m,2);p_cov(1,:)=[q1 p0];

% Density
f=zeros(m,100);
xi=zeros(m,100);
[ff,xii]=ksdensity(x_particle2);
f(1,:)=ff;xi(1,:)=xii';

% Main Loop
for i=1:m-1
  
% Truth    
 x(i+1,1)=cos(x(i,1))+sin(x(i,1))+sqrt(q1)*randn(1);
 x(i+1,2)=x(i,2)+sqrt(q2)*randn(1);
 ym(i+1)=x(i+1,1)*x(i+1,2)+sqrt(r)*randn(1);
 
 x_particle1=cos(x_particle1)+sin(x_particle1)+sqrt(q1)*randn(n,1);
 p=p+q2;
 y_cov=x_particle1.*p.*x_particle1+r;
 w_particle=w_particle.*exp(-(ym(i+1)-x_particle1.*x_particle2).^2./(2*y_cov));
 w_particle=w_particle/sum(w_particle);
 
 gain=p.*x_particle1./y_cov;
 x_particle2=x_particle2+gain.*(ym(i+1)-x_particle1.*x_particle2);
 p=(1-gain.*x_particle1).*p;
 
 xe(i+1,:)=[sum(w_particle.*x_particle1) sum(w_particle.*x_particle2)];
 x_diff_w=[w_particle.*(x_particle1-xe(i+1,1)) w_particle.*(x_particle2-xe(i+1,2))];%
 x_diff=[x_particle1-xe(i+1,1) x_particle2-xe(i+1,2)];
 p_mat=x_diff_w'*x_diff+[0 0;0 sum(w_particle.*p)];
 p_cov(i+1,:)=diag(p_mat)';
 
% Density 
 [ff,xii]=ksdensity(x_particle2);
 f(i+1,:)=ff;xi(i+1,:)=xii';

% Resample 
 [resamp_out,w_particle]=resample_pf([x_particle1 x_particle2 p],w_particle);
 x_particle1=resamp_out(:,1);
 x_particle2=resamp_out(:,2);
 p=resamp_out(:,3);
 
 
end

% Plot Densities
%surf(xii,t,f)
waterfall(xi,repmat(t,1,100),f)
%axis([-5 10 0 100 0 2])
ylabel('Time','fontsize',12)
xlabel('Sample Space','fontsize',12)
zlabel('Posterior Density','fontsize',12)

pause

% Plot Results
sig3=p_cov.^(0.5)*3;
plot(t,sig3(:,1),t,xe(:,1)-x(:,1),t,-sig3(:,1))
set(gca,'fontsize',12);
xlabel('Time')
ylabel('First State Errors')
pause
plot(t,sig3(:,2),t,xe(:,2)-x(:,2),t,-sig3(:,2))
set(gca,'fontsize',12);
xlabel('Time')
ylabel('Second State Errors')