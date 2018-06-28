function [qe,be,p_cov]=kal_attd(im,bm,wgm,r,sigu,sigv,poa,pog,av,dt,x0);
%function [qe,be,p_cov]=kal_attd(im,bm,wgm,r,sigu,sigv,poa,pog,av,dt,x0);
%
% This program determines the attitude of a satellite using a full
% Kalman filter (see Lefferts, Markley, and Shuster, JGCD Sept.-Oct. '82).
% This algorithm uses the discrete propagated error covariance.
% This version uses multiple sensors. Units are in radians and seconds.
%
% The inputs are:
%     im = inertial measurements [mx(3*s)], s = number of sensors
%     bm = body measurements  [mx(3*s)]
%    wgm = gyro measurements [mx3]
%      r = measurement covariance [(3*s)x(3*s)]
%   sigu = gyro bias noise standard deviation [3x3]
%   sigv = gyro noise standard deviation [3x3]
%    poa = initial error covariance of attitude [3x3]
%    pog = initial error covariance of gyro bias [3x3]
%     av = 1 for sensor available, 0 for not  [mxs]
%     dt = sampling interval
%     x0 = initial estimates of quaternions and biases ([q0 b0]), [1x7]
%
% The outputs are:
%     qe = estimated quaternions [mx4]
%     be = estimated gyro biases [mx3]
%  p_cov = diagonal covariances [mx6]

% Version 2, Written by John L Crassidis 8/22/94

% Constants and Conversions
x0=x0(:)';
[m,sen]=size(bm);sen=sen/3;
g=[-eye(3) zeros(3);zeros(3) eye(3)];
i500=0;

% Sensor to Body Measurements
%for mm=1:sen,
% bm(:,mm*3-2:mm*3)=(al(:,mm*3-2:mm*3)'*bm(:,mm*3-2:mm*3)')';
%end

% Pre-Allocate Space
qe=zeros(m,4);
be=zeros(m,3);
xe=zeros(m,6);
we=zeros(m,3);
p_cov=zeros(m,6);

% Initial Bias, Quaternion Estimate
be(1,:)=x0(1,5:7);
qe(1,:)=x0(1,1:4);

% Determine Initial Covariance and Discrete State Covariance 
p=[poa zeros(3);zeros(3) pog];
p_cov(1,:)=diag(p)';
qc(1:3,1:3)=(dt*sigv.^2+dt^3*sigu.^2/3);
qc(1:3,4:6)=(0.5*dt^2*sigu.^2);
qc(4:6,1:3)=qc(1:3,4:6)';
qc(4:6,4:6)=dt*sigu.^2;
qcc=g*qc*g';

% Main Loop
for i=1:m-1,

% Display When Every 500th Point Is Reached
if (i500==500), 
 disp(sprintf('      Filter has reached point %5i',i-1))
 i500=0;
end
i500=i500+1;

% Attitude Matrix
a1=qe(i,1)^2-qe(i,2)^2-qe(i,3)^2+qe(i,4)^2;
a2=2*(qe(i,1)*qe(i,2)+qe(i,3)*qe(i,4));
a3=2*(qe(i,1)*qe(i,3)-qe(i,2)*qe(i,4));
a4=2*(qe(i,1)*qe(i,2)-qe(i,3)*qe(i,4));
a5=-qe(i,1)^2+qe(i,2)^2-qe(i,3)^2+qe(i,4)^2;
a6=2*(qe(i,2)*qe(i,3)+qe(i,1)*qe(i,4));
a7=2*(qe(i,1)*qe(i,3)+qe(i,2)*qe(i,4));
a8=2*(qe(i,2)*qe(i,3)-qe(i,1)*qe(i,4));
a9=-qe(i,1)^2-qe(i,2)^2+qe(i,3)^2+qe(i,4)^2;
a=[a1 a2 a3;a4 a5 a6;a7 a8 a9];

% Output Matrix, Measurements, and Measurement Covariance
clear pbe pbe_cr h hout_tot z r1
[i1,j1]=find(av(i,:)==1);
for n=1:length(j1),
 pbe(:,n)=a*im(i,j1(n)*3-2:j1(n)*3)';
 pbe_cr=[0 -pbe(3,n) pbe(2,n);pbe(3,n) 0 -pbe(1,n);-pbe(2,n) pbe(1,n) 0];
 h(n*3-2:n*3,:)=[pbe_cr zeros(3)];
 z(n*3-2:n*3,1)=bm(i,j1(n)*3-2:j1(n)*3)';
 r1(n*3-2:n*3,n*3-2:n*3)=r(j1(n)*3-2:j1(n)*3,j1(n)*3-2:j1(n)*3);
 hout_tot(n*3-2:n*3,1)=pbe(:,n);
end

% Update State and Covariance Using Gain
k=p*h'*pinv(h*p*h'+r1);
p=(eye(6)-k*h)*p;
xe=(k*(z-hout_tot))';

% Save Output
be(i,:)=xe(1,4:6)+be(i,:);

xe(1:3)=0.5*xe(1:3);
qe11=qe(i,1)+xe(1,3).*qe(i,2)-xe(1,2).*qe(i,3)+xe(1,1).*qe(i,4);
qe22=-xe(1,3).*qe(i,1)+qe(i,2)+xe(1,1).*qe(i,3)+xe(1,2).*qe(i,4);
qe33=xe(1,2).*qe(i,1)-xe(1,1).*qe(i,2)+qe(i,3)+xe(1,3).*qe(i,4);
qe44=-xe(1,1).*qe(i,1)-xe(1,2).*qe(i,2)-xe(1,3).*qe(i,3)+qe(i,4);
qe(i,:)=[qe11 qe22 qe33 qe44]./norm([qe11 qe22 qe33 qe44]);

% Propagate Covariance
we(i,:)=wgm(i,:)-be(i,:);
w=norm(we(i,:));
wa=[0 we(i,3) -we(i,2);-we(i,3) 0 we(i,1);we(i,2) -we(i,1) 0];
phi11=eye(3)+wa*sin(w*dt)/w+wa*wa*(1-cos(w*dt))/w^2;
phi12=-(eye(3)*dt+wa*(1-cos(w*dt))/w^2+wa*wa*(w*dt-sin(w*dt))/w^3);
phi=[phi11 phi12;zeros(3) eye(3)];
p=phi*p*phi'+qcc;

% Propagate State
co=cos(0.5*w*dt);
si=sin(0.5*w*dt);
n1=we(i,1)/w;n2=we(i,2)/w;n3=we(i,3)/w;
qw1=n1*si;qw2=n2*si;qw3=n3*si;qw4=co;
om=[qw4  qw3 -qw2 qw1;-qw3  qw4  qw1 qw2;qw2 -qw1  qw4 qw3;-qw1 -qw2 -qw3 qw4];
qe(i+1,1:4)=(om*qe(i,1:4)')';
be(i+1,1:3)=be(i,1:3);


p_cov(i+1,:)=diag(p)';


% Display Error If Covariance Is Going Negative
if (length(find(p_cov(i,:)'>0))<6),
 disp(sprintf('error covariance is not positive, point %1i',i))
end

end