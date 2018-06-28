% This program compares an IMM to an MMAE appraoch to track
% a maneuvering object. 

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.6

% Written by John L. Crassidis 2/14

% Other Required Routines: none

% Time
dt=1;tf=600;t=[0:dt:tf]';m=length(t);

% Generate Truth (note: constant jerk is assumed)
x0_sim=[2000;0;10000;-15];
phi_sim=[1 dt 0 0;0 1 0 0;0 0 1 dt;0 0 0 1];
gam=[dt^2/2 0;dt 0;0 dt^2/2;0 dt]; 
u=[zeros((tf-200)/dt,2);0.075*ones((tf-400)/dt+1,2)];

x_sim=zeros(m,4);x_sim(1,:)=x0_sim(:)';
for i = 1:m-1
 x_sim(i+1,:)=(phi_sim*x_sim(i,:)'+gam*u(i,:)')';
end
x=[x_sim(:,1) x_sim(:,2) zeros(m,1) x_sim(:,3) x_sim(:,4) zeros(m,1)];

% State Transistion Matrix
phi=[1 dt dt^2/2 0 0 0
     0 1 dt 0 0 0
     0 0 1 0 0 0
     0 0 0 1 dt dt^2/2
     0 0 0 0 1 dt
     0 0 0 0 0 1];
 
% Continuous-Time Spectral Density 
g=[0 0;0 0;1 0;0 0;0 0;0 1];
q=1e-3*eye(2);

% Discrete-Time Process Noise Covariance
f=[0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 1;0 0 0 0 0 0];
big_a=[-f g*q*g';zeros(6) f']*dt;
big_b=expm(big_a);
qd=phi*big_b(1:6,7:12);

% Initial Conditions for Filters
x0=[2000;0;0;10000;-15;0];
p0=1e-12*eye(6);

% Generate Measurements
h=[1 0 0 0 0 0;0 0 0 1 0 0];
r=100^2*eye(2);

y1=x(:,1);y2=x(:,4);
ym1=y1+sqrt(r(1,1))*randn(m,1);
ym2=y2+sqrt(r(2,2))*randn(m,1);
ym=[ym1 ym2];

% IMM Filter State and Covariance Storage Variables
xe_imm=zeros(m,6);
p_imm_cov=zeros(m,6);
xe_imm1=zeros(m+1,6);xe_imm1(1,:)=x0(:)';
xe_imm2=zeros(m+1,6);xe_imm2(1,:)=x0(:)';
p_prop_imm1=p0;
p_prop_imm2=p0;

% Initial Weights for IMM
weights_imm=[0.5;0.5];weights_imm_store=zeros(m,2);weights_imm_store(1,:)=weights_imm';

% Transition Matrix for IMM
%trans=[0.95 0.05;0.05 0.95];
trans=[0.97 0.03;0.03 0.97];
%trans=[0.8 0.2;0.05 0.95];
%trans=[0.99999999999 0.00000000001;0.00000000001 0.99999999999];
%trans=eye(2);

% Initialize cbar
cbar1=trans(:,1)'*weights_imm;
cbar2=trans(:,2)'*weights_imm;
cbar=[cbar1;cbar2];

% MMAE Storage Variables
xe_mmae=zeros(m,6);
p_mmae_cov=zeros(m,6);

% MMAE Filter State and Covariance Storage Variables
xe_mmae1=zeros(m+1,6);xe_mmae1(1,:)=x0(:)';
xe_mmae2=zeros(m+1,6);xe_mmae2(1,:)=x0(:)';
p_prop_mmae1=p0;
p_prop_mmae2=p0;

% Initial Weights for MMAE
weights_mmae=[0.5;0.5];weights_mmae_store=zeros(m,2);weights_mmae_store(1,:)=weights_mmae';

% IMM Storage Variables
xe_mmae=zeros(m,6);
p_mmae_cov=zeros(m,6);

% Main Loop
for i=1:m,
 
% IMM Likelihood for Filter 1 
  e_cov_imm1=h*p_prop_imm1*h'+r;inv_e_cov_imm1=inv(e_cov_imm1);
  res_imm1=ym(i,:)'-h*xe_imm1(i,:)';
  like_imm1=1/sqrt(det(2*pi*e_cov_imm1))*exp(-0.5*res_imm1'*inv_e_cov_imm1*res_imm1);  

% IMM Likelihood for Filter 2 
  e_cov_imm2=h*p_prop_imm2*h'+r;inv_e_cov_imm2=inv(e_cov_imm2);
  res_imm2=ym(i,:)'-h*xe_imm2(i,:)';
  like_imm2=1/sqrt(det(2*pi*e_cov_imm2))*exp(-0.5*res_imm2'*inv_e_cov_imm2*res_imm2); 
   
% IMM Likelihood Vector 
  like_imm=[like_imm1;like_imm2];
  
% IMM Weights   
  weights_imm=cbar.*like_imm;
  weights_imm=weights_imm/sum(weights_imm);
  if i < m weights_imm_store(i+1,:)=weights_imm'; end
  
% IMM Filter 1 Update   
  gain_imm1=p_prop_imm1*h'*inv_e_cov_imm1;
  xe_imm1(i,:)=xe_imm1(i,:)+(gain_imm1*res_imm1)';
  p_up_imm1=(eye(6)-gain_imm1*h)*p_prop_imm1;
    
% IMM Filter 2 Update    
  gain_imm2=p_prop_imm2*h'*inv_e_cov_imm2;
  xe_imm2(i,:)=xe_imm2(i,:)+(gain_imm2*res_imm2)';
  p_up_imm2=(eye(6)-gain_imm2*h)*p_prop_imm2;
   
% IMM State and Covariance (After Measurement Update)
  xe_imm(i,:)=weights_imm(1)*xe_imm1(i,:)+weights_imm(2)*xe_imm2(i,:);
  p_imm=(p_up_imm1+(xe_imm1(i,:)-xe_imm(i,:))'*(xe_imm1(i,:)-xe_imm(i,:)))*weights_imm(1)...
       +(p_up_imm2+(xe_imm2(i,:)-xe_imm(i,:))'*(xe_imm2(i,:)-xe_imm(i,:)))*weights_imm(2);
  p_imm_cov(i,:)=diag(p_imm)';
  
% Compute Mixed Initial Conditions
  cbar1=trans(:,1)'*weights_imm;  
  if cbar1>1e-80,
   omega_11=trans(1,1)*weights_imm(1)/cbar1;
   omega_21=trans(2,1)*weights_imm(2)/cbar1;
  else
   cbar1 = 0;
   omega_11 = 0;
   omega_21 = 0;
  end
  cbar2=trans(:,2)'*weights_imm;
  if cbar2>1e-80,
   omega_12=trans(1,2)*weights_imm(1)/cbar2;
   omega_22=trans(2,2)*weights_imm(2)/cbar2;
  else
   cbar2 = 0;
   omega_12 = 0;
   omega_22 = 0;
  end
  cbar=[cbar1;cbar2];                                                     
  
% Initial Condition IMM State and Covariance
  x_imm1=omega_11*xe_imm1(i,:)+omega_21*xe_imm2(i,:);
  x_imm2=omega_12*xe_imm1(i,:)+omega_22*xe_imm2(i,:);
  p_imm1=(p_up_imm1+(xe_imm1(i,:)-x_imm1)'*(xe_imm1(i,:)-x_imm1))*omega_11...
        +(p_up_imm2+(xe_imm2(i,:)-x_imm1)'*(xe_imm2(i,:)-x_imm1))*omega_21;
  p_imm2=(p_up_imm1+(xe_imm1(i,:)-x_imm2)'*(xe_imm1(i,:)-x_imm2))*omega_12...
        +(p_up_imm2+(xe_imm2(i,:)-x_imm2)'*(xe_imm2(i,:)-x_imm2))*omega_22;
      
% IMM Filter 1 Propagation   
  xe_imm1(i+1,:)=(phi*x_imm1')';
  p_prop_imm1=phi*p_imm1*phi'+qd;   
   
% IMM Filter 2 Propagation     
  xe_imm2(i+1,:)=(phi*x_imm2')';
  p_prop_imm2=phi*p_imm2*phi;
 
% MMAE Likelihood for Filter 1 
  e_cov_mmae1=h*p_prop_mmae1*h'+r;inv_e_cov_mmae1=inv(e_cov_mmae1);
  res_mmae1=ym(i,:)'-h*xe_mmae1(i,:)';
  like_mmae1=1/sqrt(det(2*pi*e_cov_mmae1))*exp(-0.5*res_mmae1'*inv_e_cov_mmae1*res_mmae1); 
    
% MMAE Filter 1 Update   
  gain_mmae1=p_prop_mmae1*h'*inv_e_cov_mmae1;
  xe_mmae1(i,:)=xe_mmae1(i,:)+(gain_mmae1*res_mmae1)';
  p_up_mmae1=(eye(6)-gain_mmae1*h)*p_prop_mmae1;
   
% MMAE Likelihood for Filter 2 
  e_cov_mmae2=h*p_prop_mmae2*h'+r;inv_e_cov_mmae2=inv(e_cov_mmae2);
  res_mmae2=ym(i,:)'-h*xe_mmae2(i,:)';
  like_mmae2=1/sqrt(det(2*pi*e_cov_mmae2))*exp(-0.5*res_mmae2'*inv_e_cov_mmae2*res_mmae2); 

% MMAE Filter 2 Update    
  gain_mmae2=p_prop_mmae2*h'*inv_e_cov_mmae2;
  xe_mmae2(i,:)=xe_mmae2(i,:)+(gain_mmae2*res_mmae2)';
  p_up_mmae2=(eye(6)-gain_mmae2*h)*p_prop_mmae2;
   
% MMAE Likelihood Vector 
  like_mmae=[like_mmae1;like_mmae2];  
  
% Update MMAE Weights  
  weights_mmae=weights_mmae.*like_mmae;
  weights_mmae=weights_mmae/sum(weights_mmae);
  if i < m weights_mmae_store(i+1,:)=weights_mmae'; end

% MMAE Ouputs  
  xe_mmae(i,:)=weights_mmae(1)*xe_mmae1(i,:)+weights_mmae(2)*xe_mmae2(i,:);
  p_mmae=(p_up_mmae1+(xe_mmae1(i,:)-xe_mmae(i,:))'*(xe_mmae1(i,:)-xe_mmae(i,:)))*weights_mmae(1)...
        +(p_up_mmae2+(xe_mmae2(i,:)-xe_mmae(i,:))'*(xe_mmae2(i,:)-xe_mmae(i,:)))*weights_mmae(2);
  p_mmae_cov(i,:)=diag(p_mmae)';

% MMAE Filter 1 Propagation   
  xe_mmae1(i+1,:)=(phi*xe_mmae1(i,:)')';
  p_prop_mmae1=phi*p_up_mmae1*phi'+qd;       
   
% MMAE Filter 2 Propagation   
  xe_mmae2(i+1,:)=(phi*xe_mmae2(i,:)')';
  p_prop_mmae2=phi*p_up_mmae2*phi';
 
end

% Truncate State and Covariance Storage to m
xe_imm1=xe_imm1(1:m,:);xe_imm2=xe_imm2(1:m,:);

% 3-sigma Values for IMM and MMAE
sig3_imm=p_imm_cov.^(0.5)*3;
sig3_mmae=p_mmae_cov.^(0.5)*3;

%Plot Errors for First State and Weights
subplot(221)
plot(t,-sig3_mmae(:,1),t,xe_mmae(:,1)-x(:,1),t,sig3_mmae(:,1))
axis([0 600 -300 300])
set(gca,'YTick',[-300 -150 0 150 300])
set(gca,'fontsize',12)
ylabel('MMAE Errors (m)')
xlabel('Time (Sec)')

subplot(222)
plot(t,-sig3_imm(:,1),t,xe_imm(:,1)-x(:,1),t,sig3_imm(:,1))
axis([0 600 -300 300])
set(gca,'YTick',[-300 -150 0 150 300])
set(gca,'fontsize',12)
ylabel('IMM Errors (m)')
xlabel('Time (Sec)')

subplot(223)
plot(t,weights_mmae_store(:,1),t,weights_mmae_store(:,2))
set(gca,'fontsize',12)
axis([0 600 -0.5 1.5])
ylabel('MMAE Weights')
xlabel('Time (Sec)')
%legend('Filter 1','Filter2')

subplot(224)
plot(t,weights_imm_store(:,1),t,weights_imm_store(:,2))
set(gca,'fontsize',12)
axis([0 600 -0.5 1.5])
ylabel('IMM Weights')
xlabel('Time (Sec)')
%legend('Filter 1','Filter2')
