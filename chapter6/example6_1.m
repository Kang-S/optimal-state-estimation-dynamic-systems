% This example shows a simulation using a typical star camera 
% to determine the attitude of a rotating spacecraft using 
% the q method. This program provides plots of the available 
% stars and attitude errors with 3-sigma outliers.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 6.1

% Written by John L. Crassidis 10/11

% Other Required Routines: mappar.mat, starmeas.m, attm.m, i2b.m, quat_err.m, crossm.m

% Time
dt=1;tf=5400;t=[0:dt:tf]';m=length(t);

% Rotate about the r_3 vector which is the -x body axis
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
qe=zeros(m,4);pcov=zeros(m,3);

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
n=30;

b=meas(:,1:n);
bm=meas(:,31:2*n);
im=meas(:,61:3*n);
av=meas(:,91:100);
i500=0;

% Get Optimal Quaternion and Error
% Main loop
for i=1:m,

% Display When Every 500th Point is Reached
if (i500==500), 
 disp(sprintf('      Program has reached point %5i',i-1))
 i500=0;
end
i500=i500+1;

% Find Available Stars
[i1,j1]=find(av(i,:)==1);
% If Less Than Two Stars Then No Solution
if length(j1)>1, 
% Get W and V Matrices
clear w v
 for j=1:length(j1),
  w(:,j)=bm(i,j1(j)*3-2:j1(j)*3)'/norm(bm(i,j1(j)*3-2:j1(j)*3));
  v(:,j)=im(i,j1(j)*3-2:j1(j)*3)'/norm(im(i,j1(j)*3-2:j1(j)*3));
 end
 
% Form Matrices for Eigenvector Decomposition
 bigb=1/sig^2*w*v';
 s=bigb'+ bigb;
 z=[bigb(2,3)-bigb(3,2) bigb(3,1)-bigb(1,3) bigb(1,2)-bigb(2,1)];
 sigg=trace(bigb);
 k=[s-eye(3)*sigg z';z sigg];
 [ve,e]=eig(k);
 [hh,hh1]=max(diag(e));
 qe(i,:)=ve(:,hh1)';

% Covariance and 3-Sigma Outlier
ac=zeros(3);
 for j = 1:3:28,
  ac=ac-crossm(b(i,j:j+2))^2;
 end
pmat=inv(ac);
pcov(i,:)=sig^2*diag(pmat)';
 
end
end

% Quaternion Error, 3-Sigma Outliers and Total Number of Stars
qerr=quat_err(qe,q);erre=qerr(:,1:3)*2*180/pi;
sig3=pcov.^(0.5)*3*180/pi;
avm=sum(av,2);
j_less2=find(avm<2);qerr(j_less2)=0;

% Plot Available Stars
clf
stairs(t/60,avm)
set(gca,'Fontsize',12);
xlabel('Time (Min)')
ylabel('Number of Available Stars')
grid
axis([0 90 0 10])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[0 1 2 3 4 5 6 7 8 9 10])

disp(' Press any key to continue')
pause

% Plot Errors (note: tracker boresight is the z-axis)
subplot(311)
plot(t/60,sig3(:,1),'b',t/60,erre(:,1),'r',t/60,-sig3(:,1),'b');
set(gca,'Fontsize',12);
ylabel('Roll (Deg)')
grid
axis([0 90 -0.01 0.01])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.01 -0.005 0 0.005 0.01])
subplot(312)
plot(t/60,sig3(:,2),'b',t/60,erre(:,2),'r',t/60,-sig3(:,2),'b');
set(gca,'Fontsize',12);
ylabel('Pitch (Deg)')
grid
axis([0 90 -0.01 0.01])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.01 -0.005 0 0.005 0.01])
subplot(313)
plot(t/60,sig3(:,3),'b',t/60,erre(:,3),'r',t/60,-sig3(:,3),'b');
set(gca,'Fontsize',12);
ylabel('Yaw (Deg)')
grid
axis([0 90 -0.2 0.2])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.2 -0.1 0 0.1 0.2])
xlabel('Time (Min)')