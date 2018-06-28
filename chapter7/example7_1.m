% This example shows a simulation using a typical star camera 
% with gyro measurements to determine the attitude of a rotating 
% spacecraft using an extended Kalman filter. This program provides 
% plots of the available stars, the attitude errors with 3-sigma 
% outliers, and the gyro bias estimates.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 7.1

% Written by John L. Crassidis 10/11

% Other Required Routines:  mappar.mat, starmeas.m, attm.m, i2b.m, kal_attd.m, quat_err.m

% Time
dt=1;tf=5400;t=[0:dt:tf]';m=length(t);

% Rotate about the -y-body axis
w=[zeros(m,1) -inv(91.5/(2*pi)*60)*ones(m,1) zeros(m,1)];

% Initial Quaternion (note: in starmeas it is assumed that the boresight
% axis is the z axis)
% Orbit is Assumed to be in the Equatorial Plane.
% X-Body Axis in Orbit Velocity Direction, and Z-Body Axis Nadir.
% First Rotate +90 Degrees about X-Body Axis.
% Then Rotate 180 Degrees about New Z-Body Axis.
% Use a 1-3-1 Rotation with Last Angle = 0; gives q0=sqrt(2)/2)*[0;1;1;0];
q0=[cos(pi/2)*sin(-pi/4);sin(pi/2)*sin(pi/4);sin(pi/2)*cos(-pi/4);cos(pi/2)*cos(-pi/4)];
q=zeros(m,4);q(1,:)=q0(:)';

% Quaternion Over Time
for i=1:m-1,
 wnorm=norm(w(i,:));
 co=cos(0.5*wnorm*dt);
 si=sin(0.5*wnorm*dt);
 n1=w(i,1)/wnorm;n2=w(i,2)/wnorm;n3=w(i,3)/wnorm;
 qw1=n1*si;qw2=n2*si;qw3=n3*si;qw4=co;
 om=[qw4  qw3 -qw2 qw1;-qw3  qw4  qw1 qw2;qw2 -qw1  qw4 qw3;-qw1 -qw2 -qw3 qw4];
 q(i+1,1:4)=(om*q(i,1:4)')';
end

% Get Body and Inertial Vectors
% 18 arc-sec (3 sigma) or 0.005 degrees (3 sigma)
sig=87.2665/3*1e-6;
meas=starmeas(q,6,6,sig);
n=30;r=sig^2*eye(n);

b=meas(:,1:n);
bm=meas(:,31:2*n);
im=meas(:,61:3*n);
av=meas(:,91:100);

% Get True Rates (note: bias is simulated in discrete-time)
% From Reynolds, R., "Maximum Likelihood Estimation of Stability Parameters for the Standard 
% Gyroscopic Error Model," Flight Mechanics Symposium, NASA-Goddard Space Flight Center, 
% Greenbelt, MD, Oct. 2003, NASA CP-2003-212246, paper #42.
wtrue=w;
sigu=sqrt(10)*1e-10;
sigv=sqrt(10)*1e-7;
num_g=dt*[1 1];den_g=2*[1 -1];
[phi_g,gam_g,c_g,d_g]=tf2ss(num_g,den_g);
bias1=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),0.1*pi/180/3600/dt);
bias2=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),0.1*pi/180/3600/dt);
bias3=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),0.1*pi/180/3600/dt);
bias=[bias1 bias2 bias3];
wgm=wtrue+sqrt(sigv^2/dt+1/12*sigu^2*dt)*randn(m,3)+bias;

% Initial Conditions and Covariances
x0=[q(1,:)';0;0;0];
poa=(0.1*pi/180)^2*eye(3);
pog=(0.2*pi/180/3600)^2*eye(3);

% EKF (note: covariance is error quaternion)
[qe,be,pcov]=kal_attd(im,bm,wgm,r,sigu*eye(3),sigv*eye(3),poa,pog,av,dt,x0);

% Compute Attitude Errors and 3-Sigma Outlier
qerr=quat_err(qe,q);erre=qerr(:,1:3)*2*180/pi;
sig3=pcov.^(0.5)*3*180/pi;

% Plot Available Stars
clf
avm=sum(av,2);
stairs(t/60,avm)
set(gca,'Fontsize',12);
xlabel('Time (Min)')
ylabel('Number of Available Stars')
grid
axis([0 90 0 10])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[0 1 2 3 4 5 6 7 8 9 10])


disp(' Press any key to continue')
disp(' ')
pause

% Errors (note: tracker boresight is the z-axis)
subplot(311)
plot(t/60,sig3(:,1),'b',t/60,erre(:,1),'r',t/60,-sig3(:,1),'b');
set(gca,'Fontsize',12);
ylabel('Roll (Deg)')
grid
axis([0 90 -0.002 0.002])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.002 -0.001 0 0.001 0.002])
subplot(312)
plot(t/60,sig3(:,2),'b',t/60,erre(:,2),'r',t/60,-sig3(:,2),'b');
set(gca,'Fontsize',12);
ylabel('Pitch (Deg)')
grid
axis([0 90 -0.002 0.002])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.002 -0.001 0 0.001 0.002])
subplot(313)
plot(t/60,sig3(:,3),'b',t/60,erre(:,3),'r',t/60,-sig3(:,3),'b');
set(gca,'Fontsize',12);
ylabel('Yaw (Deg)')
grid
axis([0 90 -0.02 0.02])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.02 -0.01 0 0.01 0.02])
xlabel('Time (Min)')

disp(' Press any key to continue')
pause

% Biases
subplot(311)
plot(t/60,be(:,1)*180/pi*3600);
set(gca,'Fontsize',12);
ylabel('{\it x} (Deg/Hr)')
grid
axis([0 90 -0.1 0.3])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.1 0 0.1 0.2 0.3])
subplot(312)
plot(t/60,be(:,2)*180/pi*3600);
set(gca,'Fontsize',12);
ylabel('{\it y} (Deg/Hr)')
grid
axis([0 90 -0.1 0.3])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.1 0 0.1 0.2 0.3])
subplot(313)
plot(t/60,be(:,3)*180/pi*3600);
set(gca,'Fontsize',12);
ylabel('{\it z} (Deg/Hr)')
grid
axis([0 90 -0.1 0.3])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.1 0 0.1 0.2 0.3])
xlabel('Time (Min)')