% This program simulates the NED attitude and position of a vehicle.
% It uses an EKF to estimate for a moving vehicle's attitude, position 
% and velocity, as well as the gyro and accelerometer biases and 
% scale factors.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 7.2

% Written by John L. Crassidis 9/03

% Other Required Routines: attm.m, crossm.m, extract.m, kal_gps_ned_scale.m
%                          ned_fun.m, ned_fun_sim.m, quat_err.m
%                          vectnorm.m. yumadata.mat

% Time and Pre-allocation of Variables
dt=0.1;tf=480;t=[0:dt:tf]';m=length(t);x=zeros(m,10);
gravity=zeros(m,1);long=zeros(m,1);lat=zeros(m,1);height=zeros(m,1);

% Initial Conditions (Lat, Long and Height Converted to ECEF)
lat0=38*pi/180;long0=-77*pi/180;height0=0;lat(1)=lat0;long(1)=long0;height(1)=height0;

% Initial Gravity
g0=9.780327*(1+0.0053024*sin(lat0)^2-0.0000058*sin(2*lat0)^2)...
    -(3.0877e-6-0.0044e-6*sin(lat0)^2)*height0+0.0072e-12*height0^2;
gravity=zeros(m,1);gravity(1)=g0;

% Body Velocity
x_rate=[5*pi/180/60*ones((m-1)/3,1);zeros(m-(m-1)/3,1)];
y_rate=[zeros((m-1)/3,1);5*pi/180/60*ones((m-1)/3,1);zeros((m-1)/3+1,1)];
z_rate=[zeros(m-(m-1)/3,1);5*pi/180/60*ones((m-1)/3,1)];
wg=[x_rate y_rate z_rate];
%wg=zeros(m,3);

% Pre-allocation of Variables for GPS Data
number_av=zeros(m,1);
max_nit=100;
rm=zeros(m,3);r_gps=zeros(m,9);
r_gps_llh=zeros(m,9);

% GPS Data from Yuma File
load yumadata
gpsdata=yuma;

% Constants
e=0.0818;a=6378137;
mu=3.98600441e14;
w_earth=7.292115e-5;

% Measurement STD and Clock Bias in Meters
sigm=5;
bias=85000;
qbias=200;
drift=lsim(0,1,1,0,sqrt(qbias)*randn(m,1),t);

% Initial Quaternion
phi=zeros(m,1);theta=zeros(m,1);psi=zeros(m,1);
phi0=0;theta0=0;psi0=0;
phi(1)=phi0*180/pi;theta(1)=theta0*180/pi;psi(1)=psi0*180/pi;

a_ned2body0=[cos(psi0)*cos(theta0) sin(psi0)*cos(theta0) -sin(theta0)
-sin(psi0)*cos(phi0)+cos(psi0)*sin(theta0)*sin(phi0) cos(psi0)*cos(phi0)+sin(psi0)*sin(theta0)*sin(phi0) cos(theta0)*sin(phi0)
sin(psi0)*sin(phi0)+cos(psi0)*sin(theta0)*cos(phi0) -cos(psi0)*sin(phi0)+sin(psi0)*sin(theta0)*cos(phi0) cos(theta0)*cos(phi0)];

q0=extract(a_ned2body0);

% NED Acceleration
acc_n=0*ones(m,1);acc_e=0*ones(m,1);acc_d_nogravity=0.01*ones(m,1)*0;
acc_ned=zeros(m,3);
acc_ned(1,:)=[acc_n(1) acc_e(1) acc_d_nogravity(1)-g0];
acc0=a_ned2body0*acc_ned(1,:)';
acc=zeros(m,3);acc(1,:)=acc0';

% Initial Velocity
vn0=200;ve0=200;vd0=-10;

% Initial Condition
x0=[q0;lat0;long0;height0;vn0;ve0;vd0];
x=zeros(m,10);x(1,:)=x0';

% Main Simulation Loop
for i=1:m-1,
    
 f1=dt*ned_fun_sim(x(i,:),acc_ned(i,:),wg(i,:));
 f2=dt*ned_fun_sim(x(i,:)+0.5*f1',acc_ned(i,:),wg(i,:));
 f3=dt*ned_fun_sim(x(i,:)+0.5*f2',acc_ned(i,:),wg(i,:));
 f4=dt*ned_fun_sim(x(i,:)+f3',acc_ned(i,:),wg(i,:));
 x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');

 lat(i+1)=x(i+1,5);long(i+1)=x(i+1,6);height(i+1)=x(i+1,7);

 gravity(i+1)=9.780327*(1+0.0053024*sin(lat(i+1))^2-0.0000058*sin(2*lat(i+1))^2)...
              -(3.0877e-6-0.0044e-6*sin(lat(i+1))^2)*height(i+1)+0.0072e-12*height(i+1)^2;

 acc_ned(i+1,:)=[acc_n(i+1) acc_e(i+1) acc_d_nogravity(i+1)-gravity(i+1)];

 a_ned2body=attm(x(i+1,1:4));

 theta(i+1)=-asin(a_ned2body(1,3))*180/pi;
 psi(i+1)=atan2(a_ned2body(1,2),a_ned2body(1,1))*180/pi;
 phi(i+1)=atan2(a_ned2body(2,3),a_ned2body(3,3))*180/pi;

 acc(i+1,:)=(a_ned2body*acc_ned(i+1,:)')';

end

% ECEF Position
bign=a./sqrt(1-e^2*sin(lat).^2);
x_ecef=(bign+height).*cos(lat).*cos(long);
y_ecef=(bign+height).*cos(lat).*sin(long);
z_ecef=(bign*(1-e^2)+height).*sin(lat);

% Nonlinear Least Squares Solution for ECEF Position
for i=1:m,

% Matrix of Stored ECEF SV Positions
 r_ecef=zeros(28,3);

% Get Time of Applicability and Current Time
 tapp=gpsdata(4);
 tcurr=tapp+t(i);

 for k = 1:28,

% Get Orbital Data
  j=(k-1)*13+1;
  ap=gpsdata(j+8);
  ma0=gpsdata(j+9);
  inc=gpsdata(j+4);
  sma=gpsdata(j+6)^2;
  ecc=gpsdata(j+2);
  an=gpsdata(j+7);
  anr=gpsdata(j+5);

% Mean Anomaly and Greenwich Longitude
  ma=ma0+sqrt(mu/(sma^3))*(tcurr-tapp);
  om=an+anr*(tcurr-tapp)-w_earth*tcurr;
 
% Eccentric Anomaly (from "Celestial Mechanics", pg. 161) and Iteration
  big_e=ma+ecc*sin(ma)+(ecc^2/2)*sin(2*ma)+(ecc^3/24)* ...
       (9*sin(3*ma)-3*sin(ma))+(ecc^4/192)*(64 * ...
       sin(4*ma)-32*sin(2*ma));
  while abs(big_e-ecc*sin(big_e)-ma)>1e-14,
   big_e=big_e-(big_e-ecc*sin(big_e)-ma)/(1-ecc*cos(big_e));
  end
 
% True Anomaly, Argument of Latitude and Radius
  cta=(cos(big_e)-ecc)/(1-ecc*cos(big_e));
  sta=sqrt(1-ecc^2)* sin(big_e)/(1-ecc*cos(big_e));
  trueanon=atan2(sta,cta);
  arglat=trueanon+ap;
  radius=sma*(1-ecc*cos(big_e));
 
% Orbital Plane Positions
  x_prime=radius*cos(arglat);
  y_prime=radius*sin(arglat);
 
% ECEF Positions
  x_coord=x_prime*cos(om)-y_prime*cos(inc)*sin(om);
  y_coord=x_prime*sin(om)+y_prime*cos(inc)*cos(om);
  z_coord=y_prime*sin(inc);
  r_ecef(k,:)=[x_coord y_coord z_coord];
 
 end

% Vectors to SVs
 rho=r_ecef-kron(ones(28,1),[x_ecef(i) y_ecef(i) z_ecef(i)]);
 rhomag=vecnorm(rho);
 rho_n=rho./kron(rhomag,[1 1 1]);

% Find Available SVs (15 degree cutoff)
 up=[cos(lat(i))*cos(long(i)) cos(lat(i))*sin(long(i)) sin(lat(i))];
 seenang=90-acos(rho_n(:,1)*up(1)+rho_n(:,2)*up(2)+rho_n(:,3)*up(3))*180/pi;
 av=find(seenang > 15);
 available_SV=av';
 r_ecef_av=r_ecef(av,:);
 ps=rhomag(av,:);
 lm=length(ps);
 number_av(i)=lm;

% Measurements 
 psm=ps+sigm*randn(lm,1)+bias*ones(lm,1)+drift(i)*ones(lm,1);

% Nonlinear Least Squares Solution for Position
 clear re
 re(1,:)=[0 0 0 0];
 xx=re(1,:)';
 dx=1e10;
 ep_stop=1e-8;
 ii=1;

 while (norm(dx) > ep_stop), 
  rhoe=r_ecef_av-kron(ones(lm,1),[re(ii,1) re(ii,2) re(ii,3)]); 
  pe=vecnorm(rhoe)+re(ii,4);
  yres=psm-pe;
  h_nls=[-rhoe./kron(vecnorm(rhoe),[1 1 1]) ones(lm,1)];
  dx=inv(h_nls'*h_nls)*h_nls'*yres;
  xx=xx+dx;
  re(ii+1,:)=xx';
  if (ii == max_nit), break; end
  ii=ii+1;
 end

% Store ECEF Solution and Covariance
 rm(i,:)=xx(1:3)';
 pcov_nls=inv(h_nls'*h_nls)*sigm^2;
 r_gps(i,:)=[pcov_nls(1,1) pcov_nls(1,2) pcov_nls(1,3) ...
            pcov_nls(2,1) pcov_nls(2,2) pcov_nls(2,3) ... 
            pcov_nls(3,1) pcov_nls(3,2) pcov_nls(3,3)];

end

% Gyros
sigv=sqrt(10)*1e-7;sigu=sqrt(10)*1e-10;

num_g=dt*[1 1];den_g=2*[1 -1];
[phi_g,gam_g,c_g,d_g]=tf2ss(num_g,den_g);
gbias1=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),1*pi/180/3600);
gbias2=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),1*pi/180/3600);
gbias3=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),1*pi/180/3600);

gwhite1=sqrt(sigv^2/dt+1/12*sigu^2*dt)*randn(m,1);
gwhite2=sqrt(sigv^2/dt+1/12*sigu^2*dt)*randn(m,1);
gwhite3=sqrt(sigv^2/dt+1/12*sigu^2*dt)*randn(m,1);

ksg1=10e-3;ksg2=10e-3;ksg3=10e-3;

wgm1=(1+ksg1)*wg(:,1)+gbias1+gwhite1;
wgm2=(1+ksg2)*wg(:,2)+gbias2+gwhite2;
wgm3=(1+ksg3)*wg(:,3)+gbias3+gwhite3;
wgm=[wgm1 wgm2 wgm3];

% Accelerometers
sigva=9.81e-6;sigua=sqrt(3600/1e12);

num_a=dt*[1 1];den_a=2*[1 -1];
[phi_a,gam_a,c_a,d_a]=tf2ss(num_a,den_a);
abias1=dlsim(phi_a,gam_a,c_a,d_a,sigua/sqrt(dt)*randn(m,1),0.003);
abias3=dlsim(phi_a,gam_a,c_a,d_a,sigua/sqrt(dt)*randn(m,1),0.003);
abias2=dlsim(phi_a,gam_a,c_a,d_a,sigua/sqrt(dt)*randn(m,1),0.003);

awhite1=sqrt(sigva^2/dt+1/12*sigua^2*dt)*randn(m,1);
awhite2=sqrt(sigva^2/dt+1/12*sigua^2*dt)*randn(m,1);
awhite3=sqrt(sigva^2/dt+1/12*sigua^2*dt)*randn(m,1);

ksa1=5e-3;ksa2=5e-3;ksa3=5e-3;

accm1=(1+ksa1)*acc(:,1)+abias1+awhite1;
accm2=(1+ksa2)*acc(:,2)+abias2+awhite2;
accm3=(1+ksa3)*acc(:,3)+abias3+awhite3;
accm=[accm1 accm2 accm3];

% Initialize Filter
x0=[phi0 theta0 psi0 lat0 long0 height0 vn0 ve0 vd0 zeros(1,12)];
p0=zeros(21);
p0(1:3,1:3)=(1/3*pi/180)^2*eye(3);
p0(4:6,4:6)=diag([1e-6^2 1e-6^2 (20/3)^2]);
p0(7:9,7:9)=1*eye(3);
p0(10:12,10:12)=(3/3*pi/180/3600)^2*eye(3);
p0(13:15,13:15)=(0.005^2)*eye(3);
p0(16:18,16:18)=(15e-3/3)^2*eye(3);
p0(19:21,19:21)=(10e-3/3)^2*eye(3);

% Run Filters
[phi_est,theta_est,psi_est,quat_est,lat_est,long_est,height_est,vn_est,ve_est,vd_est,bg,ba,sg,sa,p_cov]=kal_gps_ins_ned_scale(rm,wgm,accm,sigu,sigv,sigua,sigva,r_gps,dt,1,p0,x0);

% Errors
sig3=p_cov.^(0.5)*3;
qerr=quat_err(quat_est,x(:,1:4));att_error=qerr(:,1:3)*2*180/pi;
bg_error=(bg-[gbias1 gbias2 gbias3])*180/pi*3600;
ba_error=ba-[abias1 abias2 abias3];
sg_error=sg-kron([ksg1 ksg2 ksg3],ones(m,1));
sa_error=sa-kron([ksa1 ksa2 ksa3],ones(m,1));

% Plot Filter Results
subplot(311)
plot(t,sig3(:,1)*180/pi,'--',t,att_error(:,1),t,-sig3(:,1)*180/pi,'b--')
axis([0 480 -1 1])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-1 -0.5 0 0.5 1]);
ylabel('Roll (deg)')
title('Kalman Filter Attitude Errors')
subplot(312)
plot(t,sig3(:,2)*180/pi,'--',t,att_error(:,2),t,-sig3(:,2)*180/pi,'b--')
axis([0 480 -1 1])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-1 -0.5 0 0.5 1]);
ylabel('Pitch (deg)')
subplot(313)
plot(t,sig3(:,3)*180/pi,'--',t,att_error(:,3),t,-sig3(:,3)*180/pi,'b--')
axis([0 480 -2 2])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-2 -1 0 1 2]);
ylabel('Yaw (deg)')
xlabel('Time (sec)')

pause

subplot(311)
plot(t,sig3(:,4)*180/pi,'--',t,(lat_est-lat)*180/pi,t,-sig3(:,4)*180/pi,'b--')
axis([0 480 -1e-4 1e-4])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-1e-4 -0.5e-4 0 0.5e-4 1e-4]);
ylabel('Latitude (deg)')
title('Kalman Filter Position Errors')
subplot(312)
plot(t,sig3(:,5)*180/pi,'--',t,(long_est-long)*180/pi,t,-sig3(:,5)*180/pi,'b--')
axis([0 480 -1e-4 1e-4])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-1e-4 -0.5e-4 0 0.5e-4 1e-4]);
ylabel('Longitude (deg)')
subplot(313)
plot(t,sig3(:,6),'--',t,height_est-height,t,-sig3(:,6),'b--')
axis([0 480 -20 20])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-20 -10 0 10 20]);
ylabel('Height (m)')
xlabel('Time (sec)')

pause

subplot(311)
plot(t,sig3(:,7),'--',t,vn_est-x(:,8),t,-sig3(:,7),'b--')
axis([0 480 -3 3])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-3 -1.5 0 1.5 3]);
ylabel('North (m/s)')
title('Kalman Filter Velocity Errors')
subplot(312)
plot(t,sig3(:,8),'--',t,ve_est-x(:,9),t,-sig3(:,8),'b--')
axis([0 480 -3 3])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-3 -1.5 0 1.5 3]);
ylabel('East (m/s)')
subplot(313)
plot(t,sig3(:,8),'--',t,vd_est-x(:,10),t,-sig3(:,8),'b--')
axis([0 480 -3 3])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-3 -1.5 0 1.5 3]);
ylabel('Down (m/s)')
xlabel('Time (sec)')

pause

subplot(311)
plot(t,sig3(:,10)*180/pi*3600,'--',t,bg_error(:,1),t,-sig3(:,10)*180/pi*3600,'b--')
axis([0 480 -5 5])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-5 -2.5 0 2.5 5]);
ylabel('x (deg/hr)')
title('Kalman Filter Gyro Bias Errors')
subplot(312)
plot(t,sig3(:,11)*180/pi*3600,'--',t,bg_error(:,2),t,-sig3(:,11)*180/pi*3600,'b--')
axis([0 480 -5 5])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-5 -2.5 0 2.5 5]);
ylabel('y (deg/hr)')
subplot(313)
plot(t,sig3(:,12)*180/pi*3600,'--',t,bg_error(:,3),t,-sig3(:,12)*180/pi*3600,'b--')
axis([0 480 -5 5])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-5 -2.5 0 2.5 5]);
ylabel('z (deg/hr)')
xlabel('Time (sec)')

pause

subplot(311)
plot(t,sig3(:,13),'--',t,ba_error(:,1),t,-sig3(:,13),'b--')
axis([0 480 -0.05 0.05])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-0.05 -0.025 0 0.025 0.05]);
ylabel('x (m/s^2)')
title('Kalman Filter Accelerometer Bias Errors')
subplot(312)
plot(t,sig3(:,14),'--',t,ba_error(:,2),t,-sig3(:,14),'b--')
axis([0 480 -0.05 0.05])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-0.05 -0.025 0 0.025 0.05]);
ylabel('y (m/s^2)')
subplot(313)
plot(t,sig3(:,15),'--',t,ba_error(:,3),t,-sig3(:,15),'b--')
axis([0 480 -0.05 0.05])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-0.05 -0.025 0 0.025 0.05]);
ylabel('z (m/s^2)')
xlabel('Time (sec)')

pause

subplot(311)
plot(t,sig3(:,16),'--',t,sg_error(:,1),t,-sig3(:,16),'b--')
axis([0 480 -0.05 0.05])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-0.05 -0.025 0 0.025 0.05]);
ylabel('x')
title('Kalman Filter Gyro Scale Factor Errors')
subplot(312)
plot(t,sig3(:,17),'--',t,sg_error(:,2),t,-sig3(:,17),'b--')
axis([0 480 -0.05 0.05])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-0.05 -0.025 0 0.025 0.05]);
ylabel('y')
subplot(313)
plot(t,sig3(:,18),'--',t,sg_error(:,3),t,-sig3(:,18),'b--')
axis([0 480 -0.05 0.05])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-0.05 -0.025 0 0.025 0.05]);
ylabel('z')
xlabel('Time (sec)')

pause

subplot(311)
plot(t,sig3(:,19),'--',t,sa_error(:,1),t,-sig3(:,19),'b--')
axis([0 480 -0.01 0.01])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-0.01 -0.005 0 0.005 0.01]);
ylabel('x')
title('Kalman Filter Accelerometer Scale Factor Errors')
subplot(312)
plot(t,sig3(:,20),'--',t,sa_error(:,2),t,-sig3(:,20),'b--')
axis([0 480 -0.01 0.01])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-0.01 -0.005 0 0.005 0.01]);
ylabel('y')
subplot(313)
plot(t,sig3(:,21),'--',t,sa_error(:,3),t,-sig3(:,21),'b--')
axis([0 480 -0.01 0.01])
set(gca,'fontsize',12)
set(gca,'Xtick',[0 60 120 180 240 300 360 420 480]);
set(gca,'Ytick',[-0.01 -0.005 0 0.005 0.01]);
ylabel('z')
xlabel('Time (sec)')