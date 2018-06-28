% This example shows a simulation involving tracking the 
% vertical position of a 747 aircraft using both the 
% alpha-beta and alpha-beta-gamma filters. It outputs 
% the 3-sigma values for each filter. This program also 
% provides plots of the position and velocity tracking 
% errors using both filters.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 7.4

% Written by John L. Crassidis 9/03

% Other Required Routines: air_fun.m

% 747 Data
cd0=0.0164;cda=0.20;cdde=0;
cy0=0;cyb=-0.9;cydr=0.120;cyda=0;
cl0=0.21;cla=4.4;clde=0.32;
cll0=0;cllb=-.160;clldr=0.008;cllda=0.013;
cm0=0;cma=-1.00;cmde=-1.30;
cn0=0;cnb=0.160;cndr=-0.100;cnda=0.0018;
cmq=-20.5;cllp=-0.340;cllr=0.130;cnp=-0.026;cnr=-0.280;

% 747 Data in SI Units
rho=.6536033;s=510.97;cbar=8.321;l=59.74;mass=2831897.6/9.81;
in=diag([24675882 44877565 67384138]);in(1,3)=1315143.1;in(3,1)=in(1,3);
g=9.81;

coef=[cd0;cda;cdde;cy0;cyb;cydr;cyda;cl0;cla;clde;cll0;cllb;clldr;cllda;cm0;cma;cmde;cn0;cnb;cndr;cnda;cmq;cllp;cllr;cnp;cnr];
other=[rho;s;cbar;l;mass;g];

% Initial Conditions
w10=0;w20=0;w30=0;
xx0=0;yy0=0;zz0=6096;
vmag=205.13;

% Trim Conditions
qtrim=0.5*rho*vmag^2;
dtrim=cla*cmde-cma*clde;
alptrim=((mass*g/qtrim/s-cl0)*cmde+cm0*clde)/dtrim;
detrim=(-cla*cm0-cma*(mass*g/qtrim/s-cl0))/dtrim;
dragtrim=(cd0+cda*alptrim+cdde*detrim)*qtrim*s;
v10=sqrt(vmag^2/(1+tan(alptrim)^2));
v20=0;
v30=v10*tan(alptrim);

% Initial Angles
phi0=0;theta0=0*alptrim;psi0=0;

% Steady-State Values
w1_ss=w10;w2_ss=w20;w3_ss=w30;
v_ss=vmag;

% True States
dt=1;
t=[0:dt:3600]';
m=length(t);
x=zeros(m,12);
x(1,1:3)=[v10 v20 v30];
x(1,4:6)=[w10 w20 w30];
x(1,7:9)=[xx0 yy0 zz0];
x(1,10:12)=[phi0 theta0 psi0];

% Control Surface Inputs and Thrust
de=detrim*ones(m,1)+1*pi/180*sin(0.01*t);
dr=0*ones(m,1);
da=0*ones(m,1);
thrust=dragtrim*ones(m,1);

% Main Loop for Aircraft Simulation
for i=1:m-1,
 f1=dt*air_fun(x(i,:),de(i),dr(i),da(i),thrust(i),coef,other,in,w1_ss,w2_ss,w3_ss,v_ss);
 f2=dt*air_fun(x(i,:)+0.5*f1',de(i),dr(i),da(i),thrust(i),coef,other,in,w1_ss,w2_ss,w3_ss,v_ss);
 f3=dt*air_fun(x(i,:)+0.5*f2',de(i),dr(i),da(i),thrust(i),coef,other,in,w1_ss,w2_ss,w3_ss,v_ss);
 f4=dt*air_fun(x(i,:)+f3',de(i),dr(i),da(i),thrust(i),coef,other,in,w1_ss,w2_ss,w3_ss,v_ss);
 x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
end

% Velocity, Angle of Attack and Sideslip
vel=x(:,1:3);
velm=(vel(:,1).^2+vel(:,2).^2+vel(:,3).^2).^(0.5);
alp=atan(vel(:,3)./vel(:,1))*180/pi;
bet=asin(vel(:,2)./velm)*180/pi;

% Angles
w=x(:,4:6)*180/pi;
pos=x(:,7:9);
phi=x(:,10)*180/pi;
theta=x(:,11)*180/pi;
psi=x(:,12)*180/pi;

% Pre-Allocate Space
pos0=pos(1,3);vel0=vel(1,3);
m=length(pos);
pose=zeros(m,1);pose(1)=pos0;
vele=zeros(m,1);vele(1)=vel0;

% Noise Parameter
sigp=10;
ym=pos(:,3)+sigp*randn(m,1);

% Gains
alp=0.1;bet=0.1;

% Alpha-Beta Filter Variables
xe_2=zeros(m,2);xe_2(1,:)=[pos0 vel0];
phi=[1 dt;0 1];h=[1 0];

% Process Noise Tuning and Covariance
q=.5;
qd=q*[1/3*dt^3 0.5*dt^2;0.5*dt^2 dt];
pcov=dare(phi',h',qd,sigp^2,zeros(2,1),eye(2));
gain=pcov*h'*inv(h*pcov*h'+sigp^2);

sig3_alp_bet=diag(pcov).^(0.5)*3
disp(' ')

% Alpha-Beta Filter
for i = 1:m-1
 xe_2(i+1,:)=(phi*xe_2(i,:)'+phi*gain*(ym(i)-xe_2(i,1)))';
end

% Alpha-Beta-Gamma Filter Variables
xe_3=zeros(m,3);xe_3(1,:)=[pos0 vel0 0];
phi=[1 dt dt^2/2;0 1 dt;0 0 1];h=[1 0 0];

% Process Noise Tuning and Covariance
q=.0001;
qd=q*[dt^5/20 dt^4/8 dt^3/6;dt^4/8 dt^3/3 dt^2/2;dt^3/6 dt^2/2 dt];
pcov=dare(phi',h',qd,sigp^2,zeros(3,1),eye(3));
gain=pcov*h'*inv(h*pcov*h'+sigp^2);

sig3_alp_bet_gam=diag(pcov).^(0.5)*3

% Alpha-Beta-Gamma Filter
for i = 1:m-1
 xe_3(i+1,:)=(phi*xe_3(i,:)'+phi*gain*(ym(i)-xe_3(i,1)))';
end

% Velocity
zvel=diff(pos(:,3))/dt;zvel(m)=zvel(m-1);

% Plot Results
subplot(221)
plot(t/60,pos(:,3)-xe_2(:,1));
set(gca,'fontsize',12)
axis([0 60 -20 20]);
set(gca,'Ytick',[-20 -10 0 10 20]);
set(gca,'Xtick',[0 20 40 60]);
ylabel('Position Error (m)')
xlabel('Time (Min)')
title('{\it \alpha}-{\it \beta} Filter')

subplot(223)
plot(t/60,zvel-xe_2(:,2));
set(gca,'fontsize',12)
axis([0 60 -3 3]);
set(gca,'Ytick',[-3 -2 -1 0 1 2 3]);
set(gca,'Xtick',[0 20 40 60]);
ylabel('Velocity Error (m/s)')
xlabel('Time (Min)')

subplot(222)
plot(t/60,pos(:,3)-xe_3(:,1));
set(gca,'fontsize',12)
axis([0 60 -20 20]);
set(gca,'Ytick',[-20 -10 0 10 20]);
set(gca,'Xtick',[0 20 40 60]);
ylabel('Position Error (m)')
xlabel('Time (Min)')
title('{\it \alpha}-{\it \beta}-{\it \gamma} Filter')

subplot(224)
plot(t/60,zvel-xe_3(:,2));
set(gca,'fontsize',12)
axis([0 60 -3 3]);
set(gca,'Ytick',[-3 -2 -1 0 1 2 3]);
set(gca,'Xtick',[0 20 40 60]);
ylabel('Velocity Error (m/s)')
xlabel('Time (Min)')
