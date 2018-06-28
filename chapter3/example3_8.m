% This example using a constrained Kalman filter to track a
% vehicle traveling down a known road with some heading angle,
% measured clockwise from due East. The state consists of the 
% north and east position, and the north and east velocity of the vehicle.
% The measurement consists of ranges to two transponders.
% For further details see the web site 
% http://academic.csuohio.edu/simond/kalmanconstrained.

% Optimal Estimation of Dynamic Systems by (2nd ed.) Crassidis and Junkins
% Example 3.8

% Modified by John L. Crassidis 10/10

% Other Required Routines: none

%function KalmanConstrained

% Time
dt=0.3;tf=300;t=[0:dt:tf]';m=length(t);

% Process Noise Covariance 
q=diag([4, 4, 1, 1]);  

% Measurement Noise Covariance 
r=diag([900, 900]); 

% Heading Angle (measured CCW from east)
theta=pi/3; 
tantheta=tan(theta);
sintheta=sin(theta);
costheta=cos(theta);

% Iniitial States
x0=[0;0;tantheta;1]*100;
x=x0;
xbar=x0; 
xe=x0; 

% Initial Estimation Error Covariance
pbar=diag([r(1,1), r(2,2), q(1,1), q(2,2)]);

% Variable used to simulate the vehicle alternately accelerating and
% decelerating, as if in traffic
acceldecelflag=1;

% Transponder Locations
rn1=0;
re1=0;
rn2=1e5*tantheta;
re2=1e5;

% System Matrix
phi = [1 0 dt 0
       0 1 0 dt
       0 0 1 0
       0 0 0 1];

% State Constraint Matrix
d = [1 -tantheta 0 0; 
     0 0 1 -tantheta];

% Initialize arrays for saving data for plotting.
x_store=zeros(m,4);x_store(1,:)=x';
xe_store=zeros(m,4);xe_store(1,:)=x';

% Covariance for Constrained Estimate
%a=eye(4)-pbar*d'*inv(d*pbar*d')*d;
%p=a*pbar*a';
%sig3=zeros(m,4);sig3(1,:)=diag(p)'.^(0.5)*3;
sig3=zeros(m,4);

% Main Loop
for i=1:m

% Measurement
 y=[(x(1)-rn1)^2 + (x(2)-re1)^2;(x(1)-rn2)^2 + (x(2)-re2)^2];
 ym=y+[sqrt(r(1,1))*randn(1);sqrt(r(2,2))*randn(1)];
      
% Constrain the Vehicle (i.e., the true state) to the Straight Road
 if abs(x(1) - tantheta * x(2)) > 0
% if abs(x(1) - tantheta * x(2)) > 2    
  x(2) = (x(2) + x(1) * tantheta) / (1 + tantheta^2);
  x(1) = x(2) * tantheta;
 end
% if abs(x(3) - tantheta * x(4)) > 0.2
 if abs(x(3) - tantheta * x(4)) > 0
  x(4) = (x(4) + x(3) * tantheta) / (1 + tantheta^2);
  x(3) = x(4) * tantheta;
 end
 x_store(i,:)=x';
 
 % Set the Known Input
 if acceldecelflag == 1
  if (x(3) > 30) | (x(4) > 30)
   acceldecelflag = -1;
  end
 else
  if (x(3) < 5) | (x(4) < 5)
   acceldecelflag = 1;
  end
 end
 u=1*acceldecelflag;

% Propagate the State
 gam=[0;0;dt*sintheta;dt*costheta];
 x=phi*x + gam*u + [sqrt(q(1,1))*randn(1);sqrt(q(2,2))*randn(1);sqrt(q(3,3))*randn(1);sqrt(q(4,4))*randn(1)];

% Unconstrained Kalman Filter Update
 h=[2*(xbar(1)-rn1) 2*(xbar(2)-re1)  0 0
    2*(xbar(1)-rn2) 2*(xbar(2)-re2) 0 0];
 gainbar=pbar*h'*inv(h*pbar*h'+r);
 pbar=(eye(4)-gainbar*h)*pbar*(eye(4)-gainbar*h)'+gainbar*r*gainbar';
 ye=[(xbar(1)-rn1)^2 + (xbar(2)-re1)^2;(xbar(1)-rn2)^2 + (xbar(2)-re2)^2];
 xbar=xbar+gainbar*(ym-ye);

% Constrained Estimate 
 gain=pbar*d'*inv(d*pbar*d');
 xe=xbar-gain*d*xbar;
 xe_store(i,:)=xe';

% Covariance for Constrained Estimate
 p=(eye(4)-gain*d)*pbar*(eye(4)-gain*d)';
 sig3(i,:)=diag(p)'.^(0.5)*3;
 
% Unconstrained Kalman Filter Propagation
 xbar=phi*xbar+gam*u;
 pbar=phi*pbar*phi'+q;

end

clf
plot(t,(d*xe_store')')
set(gca,'fontsize',12)
xlabel('Time (Sec)')
ylabel('Constraint')

pause

subplot(221)
plot(t,-sig3(:,1),t,xe_store(:,1)-x_store(:,1),t,sig3(:,1))
set(gca,'fontsize',12)
axis([0 300 -4e-4 4e-4]')
set(gca,'XTick',[0 50 100 150 200 250 300])
set(gca,'YTick',[-4e-4 -2e-4 0 2e-4 4e-4])
xlabel('Time (Sec)')
ylabel('{\it x}_1 Error and 3{\sigma} Outlier')
%ylabel('x1error')

subplot(222)
plot(t,-sig3(:,2),t,xe_store(:,2)-x_store(:,2),t,sig3(:,2))
set(gca,'fontsize',12)
axis([0 300 -2e-4 2e-4]')
set(gca,'XTick',[0 50 100 150 200 250 300])
set(gca,'YTick',[-2e-4 -1e-4 0 1e-4 2e-4])
xlabel('Time (Sec)')
ylabel('{\it x}_2 Error and 3{\sigma} Outlier')
%ylabel('x2error')

subplot(223)
plot(t,-sig3(:,3),t,xe_store(:,3)-x_store(:,3),t,sig3(:,3))
set(gca,'fontsize',12)
axis([0 300 -10 10]')
set(gca,'XTick',[0 50 100 150 200 250 300])
set(gca,'YTick',[-10 -5 0 5 10])
xlabel('Time (Sec)')
ylabel('{\it x}_3 Error and 3{\sigma} Outlier')
%ylabel('x3error')

subplot(224)
plot(t,-sig3(:,4),t,xe_store(:,4)-x_store(:,4),t,sig3(:,4))
set(gca,'fontsize',12)
axis([0 300 -5 5]')
set(gca,'XTick',[0 50 100 150 200 250 300])
set(gca,'YTick',[-5 -2.5 0 2.5 5])
xlabel('Time (Sec)')
ylabel('{\it x}_4 Error and 3{\sigma} Outlier')
%ylabel('x4error')