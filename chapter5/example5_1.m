% This example solves the single axis attitude estimation problem 
% with angle-attitude measurements and rate information from gyros
% using an RTS smoother. This program provides plots of the 
% attitude errors with 3-sigma outliers and gyro bias estimates
% from the Kalman filter and RTS smoother.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 5.1

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% True Rate
dt=1;tf=3600;t=[0:dt:tf]';m=length(t);
wtrue=0.0011;

% Gyro and Attitude Parameters
sigu=sqrt(10)*1e-10;
sigv=sqrt(10)*1e-7;
sign=17*1e-6;

% Measurements with Bias
ym=t*wtrue+sign*randn(m,1);
bias=lsim(0,1,1,0,sigu*randn(m,1),t,0.1*pi/180/3600);
wm=wtrue+sigv*randn(m,1)+bias;

% Discrete-Time Process Noise Covriance
q=[sigv^2*dt+1/3*sigu^2*dt^3 -1/2*sigu^2*dt^2;-1/2*sigu^2*dt^2 sigu^2*dt];
phi=[1 -dt;0 1];gam=[dt;0];

% Initial Covariance
poa=1e-4;
pog=1e-12;
p=[poa 0;0 pog];
pcov=zeros(m,2);pcov(1,:)=[poa pog];

% Initial Condition and H Matrix (constant)
x0=[wm(1);0];xe=zeros(m,2);xe(1,:)=x0';x=x0;
h=[1 0];

% Preallocate Variables
xup=zeros(m,2);xprop=zeros(m,2);
pup=zeros(m,3);pprop=zeros(m,3);
psmooth=zeros(m,2);
xs=zeros(m,2);

% Main Forward Loop
for i = 1:m-1

% Kalman Gain  
gain=p*h'*inv(h*p*h'+sign^2);

% Update
x=x+gain*(ym(i)-h*x);
p=[eye(2)-gain*h]*p;
xup(i+1,:)=x';
pup(i+1,:)=[p(1,1) p(1,2) p(2,2)];

% Propagate
x=phi*x+gam*wm(i);
p=phi*p*phi'+q;
xprop(i+1,:)=x';
pprop(i+1,:)=[p(1,1) p(1,2) p(2,2)];

% Store Variables
xe(i+1,:)=x';
pcov(i+1,:)=diag(p)';

end

% Initialize RTS Smoother
xs(m,:)=xe(m,:);
ps=p;

% Main Backward Loop
for i=m:-1:2,

% RTS Gain
pf_plus=[pup(i,1) pup(i,2);pup(i,2) pup(i,3)];
pf_minus=[pprop(i,1) pprop(i,2);pprop(i,2) pprop(i,3)];
gain_s=pf_plus*phi'*inv(pf_minus);

% Smoother State and Covariance
xs(i-1,:)=xup(i,:)+(gain_s*(xs(i,:)-xprop(i,:))')';
ps=pf_plus-gain_s*(pf_minus-ps)*gain_s';

% Store RTS Covariance
psmooth(i-1,:)=diag(ps)';
    
end

% 3-Sigma Outliers
sig3=pcov.^(0.5)*3;
sig3s=psmooth.^(0.5)*3;

plot(t/60,[sig3s(:,1) xs(:,1)-t*wtrue -sig3s(:,1)]*1e6)
axis([0 60 -20 20]);grid
set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');
xlabel('Time (Min)');
hh=get(gca,'Ylabel');
set(hh,'String','\fontsize{12} {Attitude Error ({\mu}rad)}');

disp(' Press any key to continue')
pause

subplot(211)
plot(t/60,xe(:,2)*180*3600/pi);grid
axis([0 60 -0.1 0.15])
set(gca,'Ytick',[-0.1 -0.05 0 0.05 0.1 0.15])
set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');
xlabel('Time (Min)');
ylabel('Bias Estimate (Deg/Hr)')

subplot(212)
plot(t/60,xs(:,2)*180*3600/pi);grid
axis([0 60 -0.1 0.15])
set(gca,'Ytick',[-0.1 -0.05 0 0.05 0.1 0.15])
set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');
xlabel('Time (Min)');
ylabel('Bias Estimate (Deg/Hr)')

% Compute Steady-State Covariances
pf_minus_steady=dare(phi',h',q,sign^2);
gainnn=pf_minus_steady*h'*inv(h*pf_minus_steady*h'+sign^2);
pf_plus_steady=[eye(2)-gainnn*h]*pf_minus_steady;

gainss=pf_plus_steady*phi'*inv(pf_minus_steady);
ps_steady=dlyap(gainss,pf_plus_steady-gainss*pf_minus_steady*gainss');