% This example solves the single axis attitude estimation problem 
% with angle-attitude measurements and rate information from gyros
% using a simple Kalman filter. This program provides plots of the 
% attitude errors with 3-sigma outliers and gyro bias estimate. It 
% also outputs the covariance solution from the Kalman filter and 
% the solution obtained by solving the discrete-time Riccati equation.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 3.3

% Written by John L. Crassidis 9/03

% Other Required Routines: none

% True Rate
dt=1;tf=3600;t=[0:dt:tf]';m=length(t);
wtrue=0.0011;

% Gyro and Attitude Parameters
sigu=sqrt(10)*1e-10;
sigv=sqrt(10)*1e-7;
sign=17*1e-6;

% Measurements with Bias (note: this is in discrete-time)
% From Reynolds, R., "Maximum Likelihood Estimation of Stability Parameters for the Standard 
% Gyroscopic Error Model," Flight Mechanics Symposium, NASA-Goddard Space Flight Center, 
% Greenbelt, MD, Oct. 2003, NASA CP-2003-212246, paper #42.
ym=t*wtrue+sign*randn(m,1);
num_g=dt*[1 1];den_g=2*[1 -1];
[phi_g,gam_g,c_g,d_g]=tf2ss(num_g,den_g);
bias=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),0.1*pi/180/3600/dt);
wm=wtrue+sqrt(sigv^2/dt+1/12*sigu^2*dt)*randn(m,1)+bias;

% Discrete-Time Process Noise Covriance
q=[sigv^2*dt+1/3*sigu^2*dt^3 -1/2*sigu^2*dt^2;-1/2*sigu^2*dt^2 sigu^2*dt];
phi=[1 -dt;0 1];gam=[dt;0];

% Initial Covariance
poa=1e-4;
pog=1e-12;
p=[poa 0;0 pog];
pcov=zeros(m,2);pcov(1,:)=[poa pog];

% Initial Condition and H Matrix (constant)
x0=[ym(1);0];xe=zeros(m,2);xe(1,:)=x0';x=x0;
h=[1 0];

% Main Loop
for i = 1:m-1

% Kalman Gain    
gain=p*h'*inv(h*p*h'+sign^2);

% Update
x=x+gain*(ym(i)-h*x);
p=[eye(2)-gain*h]*p;

% Propagate
x=phi*x+gam*wm(i);
p=phi*p*phi'+q;

% Store Variables
xe(i+1,:)=x';
pcov(i+1,:)=diag(p)';

end

% 3-Sigma Outlier
sig3=pcov.^(0.5)*3;

% Plot Results
plot(t/60,[sig3(:,1) xe(:,1)-t*wtrue -sig3(:,1)]*1e6)
axis([0 60 -20 20]);grid
set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');
xlabel('Time (Min)');
hh=get(gca,'Ylabel');
set(hh,'String','\fontsize{12} {Attitude Error ({\mu}rad)}');

disp(' Press any key to continue')
pause

plot(t/60,xe(:,2)*180*3600/pi);grid
set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');
xlabel('Time (Min)');
ylabel('Bias Estimate (Deg/Hr)')

% Show Solution with Riccati Equation
format short e
p_ric=dare(phi',h',q,sign^2,zeros(2,1),eye(2))
p