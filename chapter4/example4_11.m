% In this example a particle filter is used to estimate the posterior 
% pdf of a nonlinear discrete-time system with additive noise terms.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.11

% Written by John L. Crassidis 2/10

% Other Required Routines: resample_pf

% Time 
t=[0:1:50]';m=length(t);

% Standard Devations
sig_x=sqrt(10);sig_y=1;sig_w=sqrt(0.1);

% States and Measurements
x=zeros(m,1);x(1)=sig_x*randn(1);
ym=zeros(m,1);ym(1)=x(1)+sig_y*randn(1);

% Generate Truth
for i = 1:m-1
  x(i+1) = x(i)^2/(1+x(i)^3) + sig_w*randn(1);
  ym(i+1)= x(i+1)+ sig_y*randn(1);
end

% Particles
n=500;x_particle_opt=sig_x*randn(n,1);x_particle_bf=x_particle_opt;

% Optimal Filter
x_est_opt=zeros(m,1);
x_est_opt(1)=mean(x_particle_opt);
w_particle_opt=inv(n)*ones(n,1);

% Sigma and S Matrices are Constant
s=sig_w^2+sig_y^2;
sigma=sig_w^2-sig_w^4*inv(s);

% Density
fff=zeros(m,100);
xi=zeros(m,100);
[ff,xii]=ksdensity(x_particle_opt);
fff(1,:)=ff;xi(1,:)=xii';

% Main Loop
for i = 1:m-1
  f=x_particle_opt.^2./(1+x_particle_opt.^3);
  a=f+sigma*1/sig_y^2*(ym(i+1)-f);
  x_particle_opt=a+sqrt(sigma)*randn(n,1);
  w_particle_opt=w_particle_opt.*exp(-(ym(i+1)-f).^2/(2*s));  
  w_particle_opt=w_particle_opt/sum(w_particle_opt);
  x_est_opt(i+1)=sum(x_particle_opt.*w_particle_opt);
  
% Density  
 [ff,xii]=ksdensity(x_particle_opt);
 fff(i+1,:)=ff;xi(i+1,:)=xii';
 
end 

% Plot Results
surf(xii,t,fff)
ylabel('Time (Sec)','fontsize',12)
xlabel('Sample Space','fontsize',12)
zlabel('Posterior Density','fontsize',12)

pause

% Bootstrap Filter
x_est_bf=zeros(m,1);
x_est_bf(1)=mean(x_particle_bf);
w_particle_bf=inv(n)*ones(n,1);

% Density
fff=zeros(m,100);
xi=zeros(m,100);
[ff,xii]=ksdensity(x_particle_bf);
fff(1,:)=ff;xi(1,:)=xii';

% Main Loop
for i = 1:m-1
  x_particle_bf=x_particle_bf.^2./(1+x_particle_bf.^3)+sig_w*randn(n,1);
  w_particle_bf=w_particle_bf.*exp(-(ym(i+1)-x_particle_bf).^2/(2*sig_y^2));  
  w_particle_bf=w_particle_bf/sum(w_particle_bf);
  x_est_bf(i+1)=sum(x_particle_bf.*w_particle_bf);

 % Resample 
 [x_particle_bf,w_particle_bf]=resample_pf(x_particle_bf,w_particle_bf);

 [ff,xii]=ksdensity(x_particle_bf);
 fff(i+1,:)=ff;xi(i+1,:)=xii';
  
end 

% Plot Results
surf(xii,t,fff)
ylabel('Time (Sec)','fontsize',12)
xlabel('Sample Space','fontsize',12)
zlabel('Posterior Density','fontsize',12)