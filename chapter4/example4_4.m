% This example implements an optimal filter and two decentralized 
% filters to triagulate an object's position from range measurements. 
% The covariance intersection approach is used to fuse the estimates 
% from the decentralized filters.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.4

% Written by John L. Crassidis 2/10

% Other Required Routines: ci_fun

% Set Flag for Case (default is double_count = 0)
% double_count = 0 uses measurements 1,2 in node 1 and measurements 3,4 in node 2 
% double_count = 1 uses measurements 1,2 in node 1 and measurements 2,3,4 in node 2
double_count=0;

% Time
t=[0:0.01:10]';m=length(t);

% True Locations
x_loc=5;y_loc=5;
x1=1*ones(m,1);y1=1*t;
x2=2*t;y2=2*ones(m,1);
x3=-3*t;y3=3*ones(m,1);
%x4=3*ones(m,1);y4=-t;
x4=3*ones(m,1);y4=1*ones(m,1);
%x4=t;y4=t;

% State Transition Matrix
phi=eye(2);

% Measurement Covariance
r=diag([0.03 0.01 0.03 0.01]);

% True Outputs
ytrue1=((x1-x_loc).^2+(y1-y_loc).^2).^(0.5);
ytrue2=((x2-x_loc).^2+(y2-y_loc).^2).^(0.5);
ytrue3=((x3-x_loc).^2+(y3-y_loc).^2).^(0.5);
ytrue4=((x4-x_loc).^2+(y4-y_loc).^2).^(0.5);

% Measurements
ym1=ytrue1+sqrt(r(1,1))*randn(m,1);
ym2=ytrue2+sqrt(r(2,2))*randn(m,1);
ym3=ytrue3+sqrt(r(3,3))*randn(m,1);
ym4=ytrue4+sqrt(r(4,4))*randn(m,1);
ym=[ym1 ym2 ym3 ym4];

% Initial Condition and Covariance for Optimal Filter
x0=[4;4];p0=(2/3)^2*eye(2);p=p0;p_cov=zeros(m,2);p_cov(1,:)=diag(p)';
xe=zeros(m,2);xe(1,:)=x0';

% Initial Conditions and Covariances for Decentralized Filters
p1=p;p_cov1=zeros(m,2);p_cov1(1,:)=diag(p1)';
p2=p;p_cov2=zeros(m,2);p_cov2(1,:)=diag(p2)';
xe1=zeros(m,2);xe1(1,:)=x0';
xe2=zeros(m,2);xe2(1,:)=x0';

% Initial CI State and Covariance
p_ci=inv(0.5*inv(p1)+0.5*inv(p2));p_cov_ci=zeros(m,2);p_cov_ci(1,:)=diag(p_ci)';
xe_ci=zeros(m,2);xe_ci(1,:)=(0.5*p_ci*inv(p1)*xe1(1,:)'+0.5*p_ci*inv(p2)*xe2(1,:)')';
omega_store=zeros(m,1);omega_store(1)=0.5;

% Naive Covariance Combination
p_naive=inv(inv(p1)+inv(p2));p_cov_naive=zeros(m,2);p_cov_naive(1,:)=diag(p_naive)';

% Main Loop
for i = 1:m-1
 
% Output Estimates for Optimal Filter     
 ye1=((x1(i)-xe(i,1))^2+(y1(i)-xe(i,2))^2)^(0.5);  
 h1=-[(x1(i)-xe(i,1)) (y1(i)-xe(i,2))]/ye1^3;
 ye2=((x2(i)-xe(i,1))^2+(y2(i)-xe(i,2))^2)^(0.5);   
 h2=-[(x2(i)-xe(i,1)) (y2(i)-xe(i,2))]/ye2^3;
 ye3=((x3(i)-xe(i,1))^2+(y3(i)-xe(i,2))^2)^(0.5);   
 h3=-[(x3(i)-xe(i,1)) (y3(i)-xe(i,2))]/ye3^3;
 ye4=((x4(i)-xe(i,1))^2+(y4(i)-xe(i,2))^2)^(0.5);   
 h4=-[(x4(i)-xe(i,1)) (y4(i)-xe(i,2))]/ye4^3;
 h=[h1;h2;h3;h4];
 ye=[ye1;ye2;ye3;ye4];   

% Update for Optimal Filter 
 gain=p*h'*inv(h*p*h'+r);
 xe(i,:)=xe(i,:)+(gain*(ym(i,:)'-ye))';
 p=(eye(2)-gain*h)*p;

% Propagation for Optimal Filter 
 xe(i+1,:)=xe(i,:);
 p=phi*p*phi';
 p_cov(i+1,:)=diag(p)';
 
% Output Estimates for Decentralized Filter 1  
 ye1=((x1(i)-xe1(i,1))^2+(y1(i)-xe1(i,2))^2)^(0.5);  
 h1=-[(x1(i)-xe1(i,1)) (y1(i)-xe1(i,2))]/ye1^3;
 ye2=((x2(i)-xe1(i,1))^2+(y2(i)-xe1(i,2))^2)^(0.5);   
 h2=-[(x2(i)-xe1(i,1)) (y2(i)-xe1(i,2))]/ye2^3;
 h=[h1;h2];
 ye=[ye1;ye2];
 ym1f=ym(i,1:2);
 r1f=r(1:2,1:2);

% Update for Decentralized Filter 1  
 gain=p1*h'*inv(h*p1*h'+r1f);
 xe1(i,:)=xe1(i,:)+(gain*(ym1f'-ye))';
 p1=(eye(2)-gain*h)*p1;
 
% Propagation for Decentralized Filter 1  
 xe1(i+1,:)=xe1(i,:);
 p1=phi*p1*phi';
 p_cov1(i+1,:)=diag(p1)';
 
% Output Estimates for Decentralized Filter 2 
 ye2=((x2(i)-xe2(i,1))^2+(y2(i)-xe2(i,2))^2)^(0.5);   
 h2=-[(x2(i)-xe2(i,1)) (y2(i)-xe2(i,2))]/ye2^3;
 ye3=((x3(i)-xe2(i,1))^2+(y3(i)-xe2(i,2))^2)^(0.5);   
 h3=-[(x3(i)-xe2(i,1)) (y3(i)-xe2(i,2))]/ye3^3;
 ye4=((x4(i)-xe2(i,1))^2+(y4(i)-xe2(i,2))^2)^(0.5);   
 h4=-[(x4(i)-xe2(i,1)) (y4(i)-xe2(i,2))]/ye4^3;
 if double_count == 1
  h=[h2;h3;h4];
  ye=[ye2;ye3;ye4]; 
  ym2f=ym(i,2:4);
  r2f=r(2:4,2:4);
 else
  h=[h3;h4];
  ye=[ye3;ye4]; 
  ym2f=ym(i,3:4);
  r2f=r(3:4,3:4);
 end
 
% Update for Decentralized Filter 2   
 gain=p2*h'*inv(h*p2*h'+r2f);
 xe2(i,:)=xe2(i,:)+(gain*(ym2f'-ye))';
 p2=(eye(2)-gain*h)*p2;
    
% Propagation for Decentralized Filter 2  
 xe2(i+1,:)=xe2(i,:);
 p2=phi*p2*phi';
 p_cov2(i+1,:)=diag(p2)'; 
 
% CI Solution Using Filters 1 and 2 
 omega=fminbnd('ci_fun',0,1,[],p1,p2);
 p_ci=inv(omega*inv(p1)+(1-omega)*inv(p2));p_cov_ci(i+1,:)=diag(p_ci)';
 xe_ci(i+1,:)=(omega*p_ci*inv(p1)*xe1(i+1,:)'+(1-omega)*p_ci*inv(p2)*xe2(i+1,:)')';
 omega_store(i+1)=omega;
 
% Naive Covariance Combination
 p_naive=inv(inv(p1)+inv(p2));p_cov_naive(i+1,:)=diag(p_naive)';
 
end
    
% Get 3-Sigma Bounds
sig3=p_cov.^(0.5)*3;   
sig31=p_cov1.^(0.5)*3;
sig32=p_cov2.^(0.5)*3;
sig3_ci=p_cov_ci.^(0.5)*3; 
sig3_naive=p_cov_naive.^(0.5)*3; 

% Plot Results
clf
plot(t,omega_store)
set(gca,'fontsize',12)
set(gca,'Xtick',[0 2 4 6 8 10])
axis([0 10 -2 2])
xlabel('Time (Sec)')
ylabel('Omega')

pause

subplot(221)
plot(t,-sig31(:,1),t,xe1(:,1)-x_loc,t,sig31(:,1))
set(gca,'fontsize',12)
axis([0 10 -2 2])
set(gca,'Xtick',[0 2 4 6 8 10])
xlabel('Time (Sec)')
ylabel('Filter 1')

subplot(222)
plot(t,-sig32(:,1),t,xe2(:,1)-x_loc,t,sig32(:,1))
set(gca,'fontsize',12)
axis([0 10 -2 2])
set(gca,'Xtick',[0 2 4 6 8 10])
xlabel('Time (Sec)')
ylabel('Filter 2')

subplot(223)
plot(t,-sig3_ci(:,1),t,xe_ci(:,1)-x_loc,t,sig3_ci(:,1))
set(gca,'fontsize',12)
axis([0 10 -2 2])
set(gca,'Xtick',[0 2 4 6 8 10])
xlabel('Time (Sec)')
ylabel('Covariance Intersection')

subplot(224)
plot(t,sig3(:,1),t,sig3_naive(:,1),'--',t,sig3_ci(:,1),'-.')
set(gca,'fontsize',12)
legend('Optimal','Naive','CI')
axis([0 10 0 2])
set(gca,'Xtick',[0 2 4 6 8 10])
xlabel('Time (Sec)')
ylabel('Bounds')
