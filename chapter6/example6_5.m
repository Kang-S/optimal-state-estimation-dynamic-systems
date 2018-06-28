% The examples uses maximum likelihood estimation to identify 
% the longitudinal parameters of a simulated 747 aircraft. 
% It outputs the current iteration and parameter values. 
% This program provides plots of the measurements with the 
% best fits.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 6.5

% Written by John L. Crassidis 9/03

% Other Required Routines: aircraft_long.m, air_longfun.m

% 747 Data
cd0=0.0164;cda=0.20;cdde=0;
cl0=0.21;cla=4.4;clde=0.32;
cm0=0;cma=-1.00;cmde=-1.30;cmq=-20.5;

% Get Truth and Measurements
aircraft_long
sigalp=0.5*pi/180;
sigvel=1;
sigtheta=0.1*pi/180;
sigw=0.01*pi/180;
alpm=alp+sigalp*randn(m,1);
velm=velmag+sigvel*randn(m,1);
thetam=theta+sigtheta*randn(m,1);
wm=w+sigw*randn(m,1);

% Initial Values for Least Squares
cd0n=0.01;cd0=cd0n;
cl0n=0.1;cl0=cl0n;
cdan=0.3;cda=cdan;
clan=3;cla=clan;
cman=-0.5;cma=cman;
cm0n=0.01;cm0=cm0n;

% Step Size for Numerical Derivatives
delta=0.002;

% Set Iterations
clear xit
xit(1,:)=[cd0n cl0n cdan clan cman cm0n];


% Main Loop
for jjj = 1:10;

disp(' ')
iteration_number=jjj

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd0 Gradient
aircraft_long
alp_o=alp;
theta_o=theta;
vel_o=velmag;
w_o=w;

xcd0_p=cd0n+delta;
cd0=xcd0_p;
aircraft_long
alp_p=alp;
theta_p=theta;
vel_p=velmag;
w_p=w;

cd0=cd0n;

grad_alp=(alp_p-alp_o)/delta;
grad_theta=(theta_p-theta_o)/delta;
grad_vel=(vel_p-vel_o)/delta;
grad_w=(w_p-w_o)/delta;
gradv1=[grad_alp/sigalp;grad_theta/sigtheta;grad_vel/sigvel;grad_w/sigw];

gz1=-((alpm-alp_o)'*grad_alp/sigalp^2+...
(thetam-theta_o)'*grad_theta/sigtheta^2+...
(velm-vel_o)'*grad_vel/sigvel^2+(wm-w_o)'*grad_w/sigw^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cl0 Gradient
xcl0_p=cl0n+delta;
cl0=xcl0_p;
aircraft_long
alp_p=alp;
theta_p=theta;
vel_p=velmag;
w_p=w;

cl0=cl0n;

grad_alp=(alp_p-alp_o)/delta;
grad_theta=(theta_p-theta_o)/delta;
grad_vel=(vel_p-vel_o)/delta;
grad_w=(w_p-w_o)/delta;
gradv2=[grad_alp/sigalp;grad_theta/sigtheta;grad_vel/sigvel;grad_w/sigw];

gz2=-((alpm-alp_o)'*grad_alp/sigalp^2+...
(thetam-theta_o)'*grad_theta/sigtheta^2+...
(velm-vel_o)'*grad_vel/sigvel^2+(wm-w_o)'*grad_w/sigw^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cda Gradient
xcda_p=cdan+delta;
cda=xcda_p;
aircraft_long
alp_p=alp;
theta_p=theta;
vel_p=velmag;
w_p=w;

cda=cdan;

grad_alp=(alp_p-alp_o)/delta;
grad_theta=(theta_p-theta_o)/delta;
grad_vel=(vel_p-vel_o)/delta;
grad_w=(w_p-w_o)/delta;
gradv3=[grad_alp/sigalp;grad_theta/sigtheta;grad_vel/sigvel;grad_w/sigw];

gz3=-((alpm-alp_o)'*grad_alp/sigalp^2+...
(thetam-theta_o)'*grad_theta/sigtheta^2+...
(velm-vel_o)'*grad_vel/sigvel^2+(wm-w_o)'*grad_w/sigw^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cla Gradient
xcla_p=clan+delta;
cla=xcla_p;
aircraft_long
alp_p=alp;
theta_p=theta;
vel_p=velmag;
w_p=w;

cla=clan;

grad_alp=(alp_p-alp_o)/delta;
grad_theta=(theta_p-theta_o)/delta;
grad_vel=(vel_p-vel_o)/delta;
grad_w=(w_p-w_o)/delta;
gradv4=[grad_alp/sigalp;grad_theta/sigtheta;grad_vel/sigvel;grad_w/sigw];

gz4=-((alpm-alp_o)'*grad_alp/sigalp^2+...
(thetam-theta_o)'*grad_theta/sigtheta^2+...
(velm-vel_o)'*grad_vel/sigvel^2+(wm-w_o)'*grad_w/sigw^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cma Gradient
xcma_p=cman+delta;
cma=xcma_p;
aircraft_long
alp_p=alp;
theta_p=theta;
vel_p=velmag;
w_p=w;

cma=cman;

grad_alp=(alp_p-alp_o)/delta;
grad_theta=(theta_p-theta_o)/delta;
grad_vel=(vel_p-vel_o)/delta;
grad_w=(w_p-w_o)/delta;
gradv5=[grad_alp/sigalp;grad_theta/sigtheta;grad_vel/sigvel;grad_w/sigw];

gz5=-((alpm-alp_o)'*grad_alp/sigalp^2+...
(thetam-theta_o)'*grad_theta/sigtheta^2+...
(velm-vel_o)'*grad_vel/sigvel^2+(wm-w_o)'*grad_w/sigw^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cm0 Gradient
xcm0_p=cm0n+delta;
cm0=xcm0_p;
aircraft_long
alp_p=alp;
theta_p=theta;
vel_p=velmag;
w_p=w;

cm0=cm0n;

grad_alp=(alp_p-alp_o)/delta;
grad_theta=(theta_p-theta_o)/delta;
grad_vel=(vel_p-vel_o)/delta;
grad_w=(w_p-w_o)/delta;
gradv6=[grad_alp/sigalp;grad_theta/sigtheta;grad_vel/sigvel;grad_w/sigw];

gz6=-((alpm-alp_o)'*grad_alp/sigalp^2+...
(thetam-theta_o)'*grad_theta/sigtheta^2+...
(velm-vel_o)'*grad_vel/sigvel^2+(wm-w_o)'*grad_w/sigw^2);

% Compute Hessian
hessvec=[gradv1 gradv2 gradv3 gradv4 gradv5 gradv6];
hess=hessvec'*hessvec;
gz=[gz1;gz2;gz3;gz4;gz5;gz6];

% Update
xxn=[cd0n;cl0n;cdan;clan;cman;cm0n];
xxn=xxn-inv(hess)*gz

cd0n=xxn(1);
cl0n=xxn(2);
cdan=xxn(3);
clan=xxn(4);
cman=xxn(5);
cm0n=xxn(6);


cd0=cd0n;
cl0=cl0n;
cda=cdan;
cla=clan;
cma=cman;
cm0=cm0n;

xit(jjj+1,:)=xxn';

end

sig3_bounds=diag(inv(hess))'.^(0.5)*3

% Show Every Ten Measurements
k=[1:10:m]';

% Plot Results
subplot(221)
plot(t,alp*180/pi,t(k),alpm(k)*180/pi,'*')
axis([0 100 -2 6])
set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');  
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {{\it \alpha} (Deg)}') 
%ylabel('alpha') 
xlabel('Time (Sec)')
grid

subplot(222)
plot(t,theta*180/pi,t(k),thetam(k)*180/pi,'*')
axis([0 100 -5 10])
set(gca,'Ytick',[-5 -2.5 0 2.5 5 7.5 10]);
set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');  
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {{\it \theta} (Deg)}')  
xlabel('Time (Sec)')
grid

subplot(223)
plot(t,velmag,t(k),velm(k),'*')           
axis([0 100 180 230])
set(gca,'Ytick',[180 190 200 210 220 230])
set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');  
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {||{\bf v}|| (M/S)}') 
xlabel('Time (Sec)')
grid

subplot(224)
plot(t,w*180/pi,t(k),wm(k)*180/pi,'*')
axis([0 100 -2 4])
set(gca,'Ytick',[-2 -1 0 1 2 3 4]);
set(gca,'Fontsize',12);set(gca,'Fontname','Helvetica');  
h=get(gca,'Ylabel');
set(h,'String','\fontsize{12} {{\it \omega}_2 (Deg/Sec)}') 
xlabel('Time (Sec)')
grid