% This example uses an extended Kalman filter algorithm to 
% determine the orbit of a spacecraft from range, azimuth 
% and elevation measurements. It outputs an approximate 
% solution using the Herrick-Gibbs approach, the true values, 
% and the least squares and EKF iterations.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 7.3

% Written by John L. Crassidis 9/03

% Other Required Routines: orbitfun.m, orbit_initial.m, vecnorm.m, orbit_det.m, 
%                          orbit_kal.m, orb_prop.m

% Initial Condition and Other Variables
xi=[7000;1000;200;4;7;2];
t0=0;dt=1;tf=100;
time=[0:dt:tf]';m=length(time);
mu=398600.64;

% Get True Orbit
x=zeros(m,6);
x(1,:)=xi';
for i = 1:m-1
 f1=dt*orbitfun(x(i,:),mu);
 f2=dt*orbitfun(x(i,:)+0.5*f1',mu);
 f3=dt*orbitfun(x(i,:)+0.5*f2',mu);
 f4=dt*orbitfun(x(i,:)+f3',mu);
 x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
end

% Measurements
k=[1:10:101]';
tm=time(k);rearth=6378;
% Note phia can cause problems with atan2 and bad initial conditions
% Be very careful
phia=5*pi/180;theta0=6*pi/180;theta=theta0+2*pi/24/3600*tm;
rhou=cos(phia)*cos(theta).*(x(k,1)-rearth*cos(phia)*cos(theta)) ...
    +cos(phia)*sin(theta).*(x(k,2)-rearth*cos(phia)*sin(theta)) ...
    +sin(phia)*(x(k,3)-rearth*sin(phia));
rhoe=-sin(theta).*(x(k,1)-rearth*cos(phia)*cos(theta)) ...
     +cos(theta).*(x(k,2)-rearth*cos(phia)*sin(theta));
rhon=-sin(phia)*cos(theta).*(x(k,1)-rearth*cos(phia)*cos(theta)) ...
     -sin(phia)*sin(theta).*(x(k,2)-rearth*cos(phia)*sin(theta)) ...
     +cos(phia)*(x(k,3)-rearth*sin(phia));

% Range, Azimuth, and Elevation 
rym=vecnorm([rhou rhoe rhon]);
ym1=rym+1*randn(length(k),1);
% Note: ym2 must be from -pi to pi in this program
ym2=atan2(rhoe,rhon)+0.01*pi/180*randn(length(k),1);
ym3=asin(rhou./rym)+0.01*pi/180*randn(length(k),1);
ymcov=diag([1^2 (0.01*pi/180)^2 (0.01*pi/180)^2]);

% Initial Guess
x0=[6990;1;1;1;1;1];

% Get Solution Using Herrick-Gibbs Approach
y=[ym1(1:3) ym2(1:3) ym3(1:3)];
x2i=orbit_initial(y,tm(1:3),phia,theta(1:3));
x_herrick_gibbs=x2i'
x_true_solution=x(k(2),:)
disp(' Press any key to continue')
disp(' ')
pause

% Least Squares Solution
y=[ym1 ym2 ym3];max=30;met=1;rad_az0=phia;rad_el=theta;
[xe,xecov]=orbit_det(x0,t0,tf,dt,y,tm,rad_az0,rad_el,max,met,ymcov);

max=3;p0=1e6*eye(6);q=zeros(6);
[xe1,xecov1]=orbit_kal(x0,t0,tf,dt,y,tm,rad_az0,rad_el,max,ymcov,p0,q);

disp(' ')
xe_ls=xe

disp(' Press any key to continue')
pause

disp(' ')
xe_ekf=xe1