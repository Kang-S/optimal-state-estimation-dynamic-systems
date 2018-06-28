% This example uses a nonlinear control law combined with 
% an extended Kalman filter to maneuver a spacecraft 
% along a desired trajectory. The assumed sensors include 
% "quaternion-out" star trackers and three-axis gyros.
% This program provides plots of the roll pointing errors, 
% the angular velocity errors and the gyro drift estimates.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 8.4

% Written by John L. Crassidis 9/03

% Other Required Routines: att_fun_quat.m, crossm.m, quat_err.m, q2e.m, i2b.m

% Time and Initialize
dt=.1;t=[0:dt:6000]';m=length(t);i500=0;
x0=[sqrt(2)/2;0;0;sqrt(2)/2;0.01*0;-0.01*0;0.001*0];
x=zeros(m,7);x(1,:)=x0(:)';
xe=zeros(m,7);wm_store=zeros(m,3);qm_store=zeros(m,4);u_store=zeros(m,3);
qd_store=zeros(m,4);
inertia=[30 10 5;10 20 3;5 3 15];
inertia_inv=inv(inertia);

% Desired Quaternion
q_d=[0;0;0;1];
qd_store(1,:)=q_d';
pkiq_d=[q_d(4)*eye(3)+crossm(q_d(1:3));-q_d(1:3)'];
w_d=[0;0.0011;0];
w_ddot=[0;0;0];
q_ddot=0.5*pkiq_d*w_d;
q_dddot=0.5*pkiq_d*w_ddot-0.25*(w_d'*w_d)*q_d;

% Linear LQR Design
a=[zeros(3) eye(3);zeros(3,6)];b=[zeros(3);eye(3)];
q_weight=0.0001*eye(6);r_weight=1*eye(3);
k_gain=lqr(a,b,q_weight,r_weight);
k2=k_gain(1:3,4:6);k1=k_gain(1:3,1:3);

% Biases for Gyros
sigu=sqrt(10)*1e-10;
sigv=sqrt(10)*1e-7;
num_g=dt*[1 1];den_g=2*[1 -1];
[phi_g,gam_g,c_g,d_g]=tf2ss(num_g,den_g);
bias1=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),0.1*pi/180/3600/dt);
bias2=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),0.1*pi/180/3600/dt);
bias3=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),0.1*pi/180/3600/dt);
bias=[bias1 bias2 bias3];
qc=zeros(6);
qc(1:3,1:3)=(dt*sigv^2*eye(3)+dt^3*sigu.^2/3*eye(3));
qc(1:3,4:6)=(0.5*dt^2*sigu.^2*eye(3));
qc(4:6,1:3)=qc(1:3,4:6)';
qc(4:6,4:6)=dt*sigu.^2*eye(3);
g=[-eye(3) zeros(3);zeros(3) eye(3)];
qcc=g*qc*g';

% Noise for Star Tracker
r=(0.001*pi/180)^2*eye(3);
noise=0.5*sqrt(r(1,1))*randn(m,3);

% Kalman Filter Parameters
p_filt=[(0.001*pi/180)^2*eye(3) zeros(3);zeros(3) (0.1*pi/180/3600)^2*eye(3)];
poa=p_filt(1:3,1:3);pog=p_filt(4:6,4:6);
qe=x0(1:4);
be=zeros(3,1);
xe(1,:)=[qe' be'];

% Main Loop
for i=1:m-1,

% display when every 500th point is reached
 if (i500==500), 
  disp(sprintf('      Controller has reached point %5i',i-1))
  i500=0;
 end
 i500=i500+1;
 
% Get States
 q=x(i,1:4)';
 w=x(i,5:7)';
 
% Generate Measurements 
 wm=w+bias(i,:)'+sqrt(sigv^2/dt+1/12*sigu^2*dt)*randn(3,1);
 wm_store(i,:)=wm';
 qc=[0 -q(3) q(2);q(3) 0 -q(1);-q(2) q(1) 0];
 pkiq=[q(4)*eye(3)+qc;-q(1:3)'];
 qm=q+pkiq*noise(i,:)';qm=qm/norm(qm);
 qm_store(i,:)=qm';
 
% Update State and Covariance Using Gain
 qmm1=-qm(4)*qe(1)-qm(3)*qe(2)+qm(2)*qe(3)+qm(1)*qe(4);
 qmm2= qm(3)*qe(1)-qm(4)*qe(2)-qm(1)*qe(3)+qm(2)*qe(4);
 qmm3=-qm(2)*qe(1)+qm(1)*qe(2)-qm(4)*qe(3)+qm(3)*qe(4);
 z=2*[qmm1 qmm2 qmm3]';
 h=[eye(3) zeros(3)];
 k_filt=p_filt*h'*inv(h*p_filt*h'+r);
 p_filt=(eye(6)-k_filt*h)*p_filt;
 xee=(k_filt*z);

% Updates 
 be=xee(4:6)+be;
 xee(1:3)=0.5*xee(1:3);
 qe11=qe(1)+xee(3).*qe(2)-xee(2).*qe(3)+xee(1).*qe(4);
 qe22=-xee(3).*qe(1)+qe(2)+xee(1).*qe(3)+xee(2).*qe(4);
 qe33=xee(2).*qe(1)-xee(1).*qe(2)+qe(3)+xee(3).*qe(4);
 qe44=-xee(1).*qe(1)-xee(2).*qe(2)-xee(3).*qe(3)+qe(4);
 qe=[qe11;qe22;qe33;qe44];
 qe=qe/norm(qe); 
 we=wm-be;
  
% Controller 
 pkiq_d=[q_d(4)*eye(3)+crossm(q_d(1:3));-q_d(1:3)'];
 pkiq_ddot=[q_ddot(4)*eye(3)+crossm(q_ddot(1:3));-q_ddot(1:3)'];
 pkiq_dddot=[q_dddot(4)*eye(3)+crossm(q_dddot(1:3));-q_dddot(1:3)'];
 pkiq=[qe(4)*eye(3)+crossm(qe(1:3));-qe(1:3)'];
 om=[-crossm(we) we;-we' 0];
 
 u=crossm(we)*inertia*we + inertia*inv(pkiq_d'*pkiq) * (1/2*(we'*we)*pkiq_d' ...
  -2*pkiq_ddot'*om - k2*pkiq_d'*om - 2*pkiq_dddot' - 2*k2*pkiq_ddot' - 2*k1*pkiq_d')*qe;
 u_store(i,:)=u';
   
% Inertia Errors (if desired) 
 inertia_m= 1.0 * inertia;
 inertia_inv_m=inv(inertia_m);
 
% Call Integration Routine
 f1=dt*att_fun_quat(x(i,:),inertia_m,inertia_inv_m,u);
 f2=dt*att_fun_quat(x(i,:)+0.5*f1',inertia_m,inertia_inv_m,u);
 f3=dt*att_fun_quat(x(i,:)+0.5*f2',inertia_m,inertia_inv_m,u);
 f4=dt*att_fun_quat(x(i,:)+f3',inertia_m,inertia_inv_m,u);
 x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4'); 
   
% Propagate Covariance
 ww=norm(we);
 wa=[0 we(3) -we(2);-we(3) 0 we(1);we(2) -we(1) 0];
 phi11=eye(3)+wa*sin(ww*dt)/ww+wa*wa*(1-cos(ww*dt))/ww^2;
 phi12=-(eye(3)*dt+wa*(1-cos(ww*dt))/ww^2+wa*wa*(ww*dt-sin(ww*dt))/ww^3);
 phi=[phi11 phi12;zeros(3) eye(3)];
 p_filt=phi*p_filt*phi'+qcc; 
 
% Propagate State
 co=cos(0.5*ww*dt);
 si=sin(0.5*ww*dt);
 n1=we(1)/ww;n2=we(2)/ww;n3=we(3)/ww;
 qw1=n1*si;qw2=n2*si;qw3=n3*si;qw4=co;
 om=[qw4  qw3 -qw2 qw1;-qw3  qw4  qw1 qw2;qw2 -qw1  qw4 qw3;-qw1 -qw2 -qw3 qw4];
 qe=(om*qe);qe=qe/norm(qe);

% Store Estimated States
 xe(i+1,:)=[qe' be'];
 
% Propagate Desired Quaternion 
 ww_d=norm(w_d);
 co=cos(0.5*ww_d*dt);
 si=sin(0.5*ww_d*dt);
 n1=w_d(1)/ww_d;n2=w_d(2)/ww_d;n3=w_d(3)/ww_d;
 qw1=n1*si;qw2=n2*si;qw3=n3*si;qw4=co;
 om=[qw4  qw3 -qw2 qw1;-qw3  qw4  qw1 qw2;qw2 -qw1  qw4 qw3;-qw1 -qw2 -qw3 qw4];
 q_d=(om*q_d);q_d=q_d/norm(q_d);
 pkiq_d=[q_d(4)*eye(3)+crossm(q_d(1:3));-q_d(1:3)'];
 w_d=[0;0.0011;0];
 w_ddot=[0;0;0];
 q_ddot=0.5*pkiq_d*w_d;
 q_dddot=0.5*pkiq_d*w_ddot-0.25*(w_d'*w_d)*q_d;
 qd_store(i+1,:)=q_d';
 
end

% Get Errors
qerr=quat_err(x(:,1:4),qd_store);
[phi,theta,psi]=q2e(qerr,7);
phi=phi*180/pi;theta=theta*180/pi;psi=psi*180/pi;

% Plot Results
subplot(211)
plot(t/60,phi)
axis([0 3 -10 90])
set(gca,'Fontsize',12);
set(gca,'Xtick',[0 0.5 1 1.5 2 2.5 3])
set(gca,'Ytick',[-10 10 30 50 70 90])
ylabel('Roll (Deg)')
grid
subplot(212)
plot(t/60,phi)
axis([50 90 -0.001 0.001])
set(gca,'Fontsize',12);
set(gca,'Xtick',[50 55 60 65 70 75 80 85 90])
set(gca,'Ytick',[-1e-3 -0.5e-3 0 0.5e-3 1e-3])
ylabel('Roll (Deg)')
xlabel('Time (Min)')
grid

disp(' Press any key to continue')
disp(' ')
pause

subplot(311)
plot(t/60,xe(:,5)*180/pi*3600);
set(gca,'Fontsize',12);
ylabel('{\it x} (Deg/Hr)')
grid
axis([0 90 -15 30])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-15 0 15 30])
subplot(312)
plot(t/60,xe(:,6)*180/pi*3600);
set(gca,'Fontsize',12);
ylabel('{\it y} (Deg/Hr)')
grid
axis([0 90 -0.5 0.3])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.5 -0.3 -0.1 0.1 0.3])
subplot(313)
plot(t/60,xe(:,7)*180/pi*3600);
set(gca,'Fontsize',12);
ylabel('{\it z} (Deg/Hr)')
grid
axis([0 90 -0.1 0.4])
set(gca,'Xtick',[0 15 30 45 60 75 90])
set(gca,'Ytick',[-0.1 0 0.1 0.2 0.3 0.4])
xlabel('Time (Min)')

disp(' Press any key to continue')
disp(' ')
pause

delta_w=x(:,5:7)-i2b(qerr,kron(ones(m,1),w_d'));

subplot(311)
plot(t/60,delta_w(:,1))
axis([0 3 -0.1 0.1])
set(gca,'Fontsize',12);
set(gca,'Xtick',[0 0.5 1 1.5 2 2.5 3])
set(gca,'Ytick',[-0.1 -0.05 0 0.05 0.1])
ylabel('{\it \delta}{\it \omega}_1 (Rad/Sec)')
grid
subplot(312)
plot(t/60,delta_w(:,2))
axis([0 3 -0.0005 0.0005])
set(gca,'Fontsize',12);
set(gca,'Xtick',[0 0.5 1 1.5 2 2.5 3])
set(gca,'Ytick',[-0.0005 -0.00025 0 0.00025 0.0005])
ylabel('{\it \delta}{\it \omega}_2 (Rad/Sec)')
grid
subplot(313)
plot(t/60,delta_w(:,3))
axis([0 3 -0.001 0.002])
set(gca,'Fontsize',12);
set(gca,'Xtick',[0 0.5 1 1.5 2 2.5 3])
set(gca,'Ytick',[-0.001 0 0.001 0.002])
ylabel('{\it \delta}{\it \omega}_3 (Rad/Sec)')
grid
xlabel('Time (Min)')