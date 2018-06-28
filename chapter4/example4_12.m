% In this example a bootstrap filter is used to estimate the posterior 
% pdf of a nonlinear discrete-time system with additive noise terms.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.12

% Written by John L. Crassidis 2/10

% Other Required Routines: resample_pf

% Time 
m=51;t=[0:1:50]';

% Standard Devations
sig_x=sqrt(5);sig_y=1;sig_w=sqrt(10);

% States and Measurements
x=zeros(m,1);x(1)=sig_x*randn(1);x_est=zeros(m,1);
ym=zeros(m,1);ym(1)=x(1)^2/20+sig_y*randn(1);

% Particles
n=500;x_particle=sig_x*randn(n,1);
x_est(1)=mean(x_particle);
w_particle=inv(n)*ones(n,1);

% Density
f=zeros(m,100);
xi=zeros(m,100);
[ff,xii]=ksdensity(x_particle);
f(1,:)=ff;xi(1,:)=xii';

% Main Loop
for i = 1:m-1
  x(i+1) = x(i)/2 + 25*x(i)/(1+x(i)^2) + 8*cos(1.2*i) + sig_w*randn(1);
  ym(i+1)= x(i+1)^2/20+ sig_y*randn(1);
  x_particle = x_particle/2 +25*x_particle./(1+x_particle.^2) + ...
                      8*cos(1.2*i) + sig_w*randn(n,1);
  w_particle = w_particle.*exp(-(ym(i+1)-x_particle.^2/20).^2/(2*sig_y^2));  
  w_particle = w_particle/sum(w_particle);
  x_est(i+1)=sum(x_particle.*w_particle);

 % Resample
  [x_particle,w_particle]=resample_pf(x_particle,w_particle);
  
 % Density 
  [ff,xii]=ksdensity(x_particle);
  f(i+1,:)=ff;xi(i+1,:)=xii';
  
end 

% Plot Results
%surf(xii,t,f)
waterfall(xi,repmat(t,1,100),f)
ylabel('Time','fontsize',12)
xlabel('Sample Space','fontsize',12)
zlabel('Posterior Density','fontsize',12)
