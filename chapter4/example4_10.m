% In this example a GSF is used to estimate the posterior pdf of
% a nonlinear discrete-time system with additive noise terms.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.10

% Written by John L. Crassidis 2/10

% Other Required Routines: none

% Time 
t=[0:1:50]';m=length(t);

% Standard Devations
r=1;q=10;

% States and Measurements
x=zeros(m,1);x(1)=0;
ym=zeros(m,1);ym(1)=x(1)^2/20+sqrt(r)*randn(1);

n=501;
weights=1/n;
xe_filters=linspace(-5,5,n)';
p=linspace(0.1,5,n)';

%n=1;xe_filters=0;p=0.1;

% Pre-allocate Space
xe=zeros(m,1);xe(1)=sum(xe_filters)/n;
p_cov=zeros(m,1);p_cov(1)=sum(p)/n;

% Main Loop
for i = 1:m

% Truth and Measurements    
  if i < m  
   x(i+1) = x(i)/2 + 25*x(i)/(1+x(i)^2) + 8*cos(1.2*i) + sqrt(q)*randn(1);
   ym(i+1)=x(i+1)^2/20+sqrt(r)*randn(1);
  end
 
 % Filters and Weights 
  h=xe_filters/10;
  bige=h.^2.*p+r;
  gain=p.*h.*bige.^(-1);
  e_res=ym(i)-xe_filters.^2/20;
  weights=weights./((2*pi*bige).^(0.5)).*exp(-0.5*e_res.*e_res./bige);
  weights=weights/sum(weights);
  
  xe_filters=xe_filters+gain.*e_res;
  p=(1-gain.*h).*p;
  
 % Estimates 
  xe(i)=sum(weights.*xe_filters);
  p_cov(i)=sum(weights.*(xe_filters-xe(i)).^2+p);
  
  phi=0.5+25./(1+xe_filters.^2)-50*xe_filters.^2./(1+xe_filters.^2).^2;
  xe_filters=xe_filters/2+25*xe_filters./(1+xe_filters.^2)+8*cos(1.2*i);
  p=phi.*p.*phi+q;
   
end 

% Plot Results
sig3=p_cov.^(0.5);

plot(t,sig3,t,x-xe,t,-sig3)
set(gca,'fontsize',12)
xlabel('Time')
ylabel('State Errors')
