function [phi_est,theta_est,psi_est,quat_est,lat_est,long_est,height_est,vn_est,ve_est,vd_est,bg,ba,sg,sa,p_cov]=kal_gps_ins_ned_scale(rm,wgm,accm,sigu,sigv,sigua,sigva,r_gps,dt,gps_up,p0,x0);
% [phi_est,theta_est,psi_est,quat_est,lat_est,long_est,height_est,vn_est,ve_est,vd_est,bg,ba,sg,sa,p_cov]=kal_gps_ins_ned_scale(rm,wgm,accm,sigu,sigv,sigua,sigva,r_gps,dt,gps_up,p0,x0);
%
% This program determines the attitude, position and velocity of a vehicle
% using a (loosely) combined GPS/INS with an extended Kalman filter.
% This algorithm uses the discrete propagated error covariance.
% Units are in radians, meters and seconds.
%
% The inputs are:
%      rm = ECEF measurements [mx3] (m)
%     wgm = gyro measurements [mx3]
%    accm = accelerometer measurements [mx3]
%    sigu = gyro bias noise standard deviation [1x1]
%    sigv = gyro noise standard deviation [1x1]
%   sigua = accelerometer bias noise standard deviation [1x1]
%   sigva = accelerometer noise standard deviation [1x1]
%   r_gps = measurement covariance at each time for llh (rcov = [r(:,1) r(:,2) r(:,3)]) [mx9]
%  gps_up = GPS update = dt_gps/dt_other_sensors [1x1]
%      dt = sampling interval
%      p0 = initial covariance (attitude error in rad^2, NED pos, NED vel, gyro bias, acc bias, scale factors) [21x21]
%      x0 = initial estimates of attitude in deg (roll, pitch, yaw), NED postion, NED velocity, biases, scale [21x1]
%
% The outputs are:
%    phi_est = estimated roll angle wrt NED (deg) [mx1]
%  theta_est = estimated pitch angle wrt NED (deg) [mx1]
%    psi_est = estimated yaw angle wrt NED (deg) [mx1]
%   quat_est = quaternion estimate wrt NED [mx4]
%    lat_est = estimated latitude (rad) [mx1]
%   long_est = estimated longitude (rad) [mx1]
% height_est = estimated height (meters) [mx1]
%     vn_est = estimated North velocity (m/s) [mx1]
%     ve_est = estimated East velocity (m/s) [mx1]
%     vd_est = estimated Down velocity (m/s) [mx1]
%         bg = estimated gyro biases (rad/sec) [mx3]
%         ba = estimated accelerometer biases (m/s^2) [mx3]
%         sg = estimated gyro scale factors [mx3]
%         sa = estimated accelerometer scale factors [mx3]
%      p_cov = diagonal covariances [mx21]


% Variables and Initial Conditions
m=length(wgm);
x0=x0(:);
phi_est=zeros(m,1);theta_est=zeros(m,1);psi_est=zeros(m,1);
lat_est=zeros(m,1);long_est=zeros(m,1);height_est=zeros(m,1);
vn_est=zeros(m,1);ve_est=zeros(m,1);vd_est=zeros(m,1);
bg=zeros(m,3);ba=zeros(m,3);sg=zeros(m,3);sa=zeros(m,3);
phi0=x0(1)*pi/180;theta0=x0(2)*pi/180;psi0=x0(3)*pi/180;
phi_est(1)=x0(1);theta_est(1)=x0(2);psi_est(1)=x0(3);
lat_est(1)=x0(4);long_est(1)=x0(5);height_est(1)=x0(6);
vn_est(1)=x0(7);ve_est(1)=x0(8);vd_est(1)=x0(9);
bg(1,:)=x0(10:12)';ba(1,:)=x0(13:15)';sg(1,:)=x0(16:18)';sa(1,:)=x0(19:21)';
quat_est=zeros(m,4);

% Initial Quaternion
a_ned2body0=[cos(psi0)*cos(theta0) sin(psi0)*cos(theta0) -sin(theta0)
-sin(psi0)*cos(phi0)+cos(psi0)*sin(theta0)*sin(phi0) cos(psi0)*cos(phi0)+sin(psi0)*sin(theta0)*sin(phi0) cos(theta0)*sin(phi0)
sin(psi0)*sin(phi0)+cos(psi0)*sin(theta0)*cos(phi0) -cos(psi0)*sin(phi0)+sin(psi0)*sin(theta0)*cos(phi0) cos(theta0)*cos(phi0)];
qe=extract(a_ned2body0);quat_est(1,:)=qe';

% Covariance
p_cov=zeros(m,21);
qcov=[sigv^2*eye(3) zeros(3) zeros(3) zeros(3)
    zeros(3) sigu^2*eye(3) zeros(3) zeros(3)
    zeros(3) zeros(3) sigva^2*eye(3) zeros(3)
    zeros(3) zeros(3) zeros(3) sigua^2*eye(3)];
p=p0;p_cov(1,:)=diag(p0)';

% Constants
a=6378137;
e=0.0818;
we=7.292155e-5;

k_up=0;
i500=0;

% Main Loop
for i=1:m-1

% Display When Every 500th Point is Reached
if (i500==500), 
 disp(sprintf('      Program has reached point %5i',i-1))
 i500=0;
end
i500=i500+1;    

% Update
if k_up == 0,
 
cos_lat=cos(lat_est(i));cos_long=cos(long_est(i));
sin_lat=sin(lat_est(i));sin_long=sin(long_est(i));
bign=a/sqrt(1-e^2*sin_lat^2);
xe_est=(bign+height_est(i))*cos_lat*cos_long;
ye_est=(bign+height_est(i))*cos_lat*sin_long;
ze_est=(bign*(1-e^2)+height_est(i))*sin_lat;
dn_dlat=a*e^2*sin_lat*cos_lat/(1-e^2*sin_lat^2)^1.5;
h11=dn_dlat*cos_lat*cos_long-(bign+height_est(i))*sin_lat*cos_long;
h12=-(bign+height_est(i))*cos_lat*sin_long;
h13=cos_lat*cos_long;
h21=dn_dlat*cos_lat*sin_long-(bign+height_est(i))*sin_lat*sin_long;
h22=(bign+height_est(i))*cos_lat*cos_long;
h23=cos_lat*sin_long;
h31=dn_dlat*(1-e^2)*sin_lat+(bign*(1-e^2)+height_est(i))*cos_lat;
h33=sin_lat;
h_sen=[h11 h12 h13;h21 h22 h23;h31 0 h33];
h=[zeros(3) h_sen zeros(3,15)];

rcov=[r_gps(i,1) r_gps(i,2) r_gps(i,3);r_gps(i,4) r_gps(i,5) r_gps(i,6);r_gps(i,7) r_gps(i,8) r_gps(i,9)];
gain=p*h'*inv(h*p*h'+rcov);
xup=(gain*(rm(i,:)'-[xe_est;ye_est;ze_est]))';

%p=(eye(21)-gain*h)*p;
p=(eye(21)-gain*h)*p*(eye(21)-gain*h)'+gain*rcov*gain';

omda=[-crossm(xup(1:3)) xup(1:3)';-xup(1:3) 0];
qe=qe+(0.5*omda*qe);
qe=qe/norm(qe);
lat_est(i)=lat_est(i)+xup(4);
long_est(i)=long_est(i)+xup(5);
height_est(i)=height_est(i)+xup(6);
vn_est(i)=vn_est(i)+xup(7);
ve_est(i)=ve_est(i)+xup(8);
vd_est(i)=vd_est(i)+xup(9);
bg(i,:)=bg(i,:)+xup(10:12);
ba(i,:)=ba(i,:)+xup(13:15);
sg(i,:)=sg(i,:)+xup(16:18);
sa(i,:)=sa(i,:)+xup(19:21);

end

k_up=k_up+1; if k_up==gps_up, k_up=0; end 

% Define Variables (need to define again because of the update)
sin_lat=sin(lat_est(i));
cos_lat=cos(lat_est(i));

r_lat=a*(1-e^2)/((1-e^2*sin_lat^2)^(1.5));
r_long=a/((1-e^2*sin_lat^2)^(0.5));

den_lat=r_lat+height_est(i);den_lat2=den_lat^2;
den_long=r_long+height_est(i);den_long2=den_long^2;

dr_lat=3*a*(1-e^2)*e^2*sin_lat*cos_lat/(1-e^2*sin_lat^2)^(5/2);
dr_long=a*e^2*sin_lat*cos_lat/(1-e^2*sin_lat^2)^(3/2);

% Attitude Matrix
a_ned2body=attm(qe);

% Estimate Rate and Acc.
wge=(eye(3)-diag(sg(i,:)))*(wgm(i,:)'-bg(i,:)');
acce=(eye(3)-diag(sa(i,:)))*(accm(i,:)'-ba(i,:)');

% Gravity 
g=9.780327*(1+0.0053024*sin(lat_est(i))^2-0.0000058*sin(2*lat_est(i))^2)...
    -(3.0877e-6-0.0044e-6*sin_lat^2)*height_est(i)+0.0072e-12*height_est(i)^2;

% Partials
f11=-crossm(wge);

domega_dp=[-we*sin_lat-ve_est(i)/den_long2*dr_long 0 -ve_est(i)/den_long2;
  vn_est(i)/den_lat2*dr_lat 0 vn_est(i)/den_lat2;
 -we*cos_lat-ve_est(i)/den_long/cos_lat^2+ve_est(i)*sin_lat/den_long2/cos_lat*dr_long 0 ve_est(i)*sin_lat/den_long2/cos_lat];
f12=-a_ned2body*domega_dp;

domega_dv=[0 1/den_long 0;-1/den_lat 0 0;0 -sin_lat/den_long/cos_lat 0];
f13=-a_ned2body*domega_dv;

f14=-(eye(3)-diag(sg(i,:)));

f16=-diag(wgm(i,:)'-bg(i,:)');

f22=[-vn_est(i)/den_lat2*dr_lat 0 -vn_est(i)/den_lat^2;
    -ve_est(i)/den_long2/cos_lat*dr_long+ve_est(i)*sin_lat/den_long/cos_lat^2 0 -ve_est(i)/den_long2/cos_lat
    0 0 0];
f23=[1/den_lat 0 0;0 1/den_long/cos_lat 0;0 0 -1];

f31=-a_ned2body'*crossm(acce);

dg_dlat=9.780327*[0.0106048*sin_lat*cos_lat-0.0000464*(sin_lat*cos_lat^3-sin_lat^3*cos_lat)]+0.0088e-6*height_est(i)*sin_lat*cos_lat;
dg_dheight=-3.0877e-6+0.0044e-6*sin_lat^2+0.0144e-12*height_est(i);


y11=-ve_est(i)^2/den_long/cos_lat^2+ve_est(i)^2*sin_lat/den_long2/cos_lat*dr_long-2*we*ve_est(i)*cos_lat-vn_est(i)*vd_est(i)/den_lat2*dr_lat;
y13=ve_est(i)^2*sin_lat/den_long2/cos_lat-vn_est(i)*vd_est(i)/den_lat2;
y21=ve_est(i)*vn_est(i)/den_long/cos_lat^2-ve_est(i)*vn_est(i)*sin_lat/den_long2/cos_lat*dr_long+2*we*(vn_est(i)*cos_lat-vd_est(i)*sin_lat)-ve_est(i)*vd_est(i)/den_long2*dr_long;
y23=-ve_est(i)*vn_est(i)*sin_lat/den_long2/cos_lat-ve_est(i)*vd_est(i)/den_long2;
y31=ve_est(i)^2/den_long2*dr_long+vn_est(i)^2/den_lat2*dr_lat+2*we*ve_est(i)*sin_lat+dg_dlat;
y33=ve_est(i)^2/den_long2+vn_est(i)^2/den_lat2+dg_dheight;
f32=[y11 0 y13;y21 0 y23;y31 0 y33];

z11=vd_est(i)/den_lat;
z12=-2*ve_est(i)*sin_lat/den_long/cos_lat-2*we*sin_lat;
z13=vn_est(i)/den_lat;
z21=ve_est(i)*sin_lat/den_long/cos_lat+2*we*sin_lat;
z22=vn_est(i)*sin_lat/den_long/cos_lat+vd_est(i)/den_long;
z23=ve_est(i)/den_long+2*we*cos_lat;
z31=-2*vn_est(i)/den_lat;
z32=-2*ve_est(i)/den_long-2*we*cos_lat;
f33=[z11 z12 z13;z21 z22 z23;z31 z32 0];

f35=-a_ned2body'*(eye(3)-diag(sa(i,:)));

f37=-a_ned2body'*diag(accm(i,:)'-ba(i,:)');


% Covariance Propagation
g_mat=[-(eye(3)-diag(sg(i,:))) zeros(3,9);zeros(3,12);zeros(3,6) -a_ned2body'*(eye(3)-diag(sa(i,:))) zeros(3);zeros(3) eye(3) zeros(3,6);zeros(3,9) eye(3);zeros(6,12)];

f_mat=[f11 f12 f13 f14 zeros(3) f16 zeros(3)
    zeros(3) f22 f23 zeros(3,12)
    f31 f32 f33 zeros(3) f35 zeros(3) f37
    zeros(12,21)];

biga=dt*[-f_mat g_mat*qcov*g_mat';zeros(21) f_mat'];
bigb=expm(biga);
phi_state=bigb(22:42,22:42)';
qd_cov=phi_state*bigb(1:21,22:42);

p=phi_state*p*phi_state'+qd_cov;

%pvec=p(:)';
%fp1=dt*covfun(pvec,f_mat,g_mat*qcov*g_mat');
%fp2=dt*covfun(pvec+0.5*fp1',f_mat,g_mat*qcov*g_mat');
%fp3=dt*covfun(pvec+0.5*fp2',f_mat,g_mat*qcov*g_mat');
%fp4=dt*covfun(pvec+fp3',f_mat,g_mat*qcov*g_mat');
%pvec=pvec+1/6*(fp1'+2*fp2'+2*fp3'+fp4');
%p(:)=pvec;

% State Propagation
x_state=[qe' lat_est(i) long_est(i) height_est(i) vn_est(i) ve_est(i) vd_est(i)];

f1=dt*ned_fun(x_state,acce,wge);
f2=dt*ned_fun(x_state+0.5*f1',acce,wge);
f3=dt*ned_fun(x_state+0.5*f2',acce,wge);
f4=dt*ned_fun(x_state+f3',acce,wge);
x_state=x_state+1/6*(f1'+2*f2'+2*f3'+f4');

qe=x_state(1:4)';
lat_est(i+1)=x_state(5);
long_est(i+1)=x_state(6);
height_est(i+1)=x_state(7);
vn_est(i+1)=x_state(8);
ve_est(i+1)=x_state(9);
vd_est(i+1)=x_state(10);
bg(i+1,:)=bg(i,:);
ba(i+1,:)=ba(i,:);
sg(i+1,:)=sg(i,:);
sa(i+1,:)=sa(i,:);

% Euler Angles
a_ned2body_new=attm(qe);
theta_est(i+1)=-asin(a_ned2body_new(1,3))*180/pi;
psi_est(i+1)=atan2(a_ned2body_new(1,2),a_ned2body_new(1,1))*180/pi;
phi_est(i+1)=atan2(a_ned2body_new(2,3),a_ned2body_new(3,3))*180/pi;
quat_est(i+1,:)=qe';

p_cov(i+1,:)=diag(p)';


end