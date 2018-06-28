% The examples uses an extended Kalman filter to identify 
% the longitudinal parameters of a simulated 747 aircraft. 
% This program provides plots of the parameter estimates 
% and pitch error with 3-sigma outliers.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 7.5

% Written by John L. Crassidis 9/03

% Other Required Routines: air_longfuns.m, air_covfuns.m, air_longfuns_ekf.m

% 747 Data
cd0=0.0164;cda=0.20;cdde=0;
cl0=0.21;cla=4.4;clde=0.32;
cm0=0;cma=-1.00;cmde=-1.30;cmq=-20.5;

% 747 Data, Used for Trim Values
cd00=0.0164;cda0=0.20;cdde0=0;
cl00=0.21;cla0=4.4;clde0=0.32;
cm00=0;cma0=-1.00;cmde0=-1.30;cmq0=-20.5;

% 747 Data in SI Units
rho=.6536033;s=510.97;cbar=8.321;l=59.74;mass=2831897.6/9.81;
in=diag([24675882 44877565 67384138]);in(1,3)=1315143.1;in(3,1)=in(1,3);
g=9.81;

coef=[cd0;cda;cdde;cl0;cla;clde;cm0;cma;cmde;cmq];
other=[rho;s;cbar;l;mass;g];

% Initial Conditions
v10=205.13;v30=0;
w20=0;
xx0=0;zz0=6096;

% Trim Conditions
qtrim=0.5*rho*norm([v10 v30])^2;
dtrim=cla0*cmde0-cma0*clde0;
alptrim=((mass*g/qtrim/s-cl00)*cmde0+cm00*clde0)/dtrim;
detrim=(-cla0*cm00-cma0*(mass*g/qtrim/s-cl00))/dtrim;
dragtrim=(cd00+cda0*alptrim+cdde0*detrim)*qtrim*s;

% Initial Angles
theta0=0;

% Steady-State Values
w2_ss=w20;
v_ss=sqrt(v10^2+v30^2);

% dt and EKF Update Time (can be changed)
dt=0.01;
update_sec=0.1;
kk=1;kup=update_sec/dt;
t=[0:dt:30]';

% True States
m=length(t);
x=zeros(m,4);
x(1,:)=[v10 v30 w20 theta0];
i500=0;

% Control Surface Inputs and Thrust
thrust=dragtrim*ones(m,1);
de=detrim*ones(m,1);de(1:5)=detrim-ones(5,1)*1*pi/180;

% Main Loop for Aircraft Simulation
for i=1:m-1,
 f1=dt*air_longfuns(x(i,:),de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
 f2=dt*air_longfuns(x(i,:)+0.5*f1',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
 f3=dt*air_longfuns(x(i,:)+0.25*f1'+0.25*f2',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
 f4=dt*air_longfuns(x(i,:)-f2'+2*f3',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
 f5=dt*air_longfuns(x(i,:)+7/27*f1'+10/27*f2'+1/27*f4',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
 f6=dt*air_longfuns(x(i,:)+28/625*f1'-1/5*f2'+546/625*f3'+54/625*f4'-378/625*f5',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
 x(i+1,:)=x(i,:)+(1/24*f1'+5/48*f4'+27/56*f5'+125/336*f6');
end

% Velocity, Angle of Attack and Pitch
vel=x(:,1:2);
velmag=(vel(:,1).^2+vel(:,2).^2).^(0.5);
alp=atan(vel(:,2)./vel(:,1));
w=x(:,3);
theta=x(:,4);

% Measurements
sigalp=0.5*pi/180;
sigvel=1;
sigtheta=0.1*pi/180;
sigw=0.01*pi/180;
alpm=alp+sigalp*randn(m,1);
velm=velmag+sigvel*randn(m,1);
thetam=theta+sigtheta*randn(m,1);
wm=w+sigw*randn(m,1);

ym=[alpm velm wm thetam];
rcov=diag([sigalp^2 sigvel^2 sigw^2 sigtheta^2]);

% Initialize Variables
xe=zeros(m,7);xe(1,1:4)=x(1,1:4);xe(1,5)=0.01;xe(1,6)=0.1;xe(1,7)=0.01;
pcov=zeros(7);pcov(1,1)=1;pcov(2,2)=1;pcov(3,3)=0.01;pcov(4,4)=0.01;
pcov(5,5)=0.1^2;pcov(6,6)=1^2;pcov(7,7)=0.1^2;
pcov=0.00001*eye(7);pcov(5,5)=1^2;pcov(6,6)=1^2;pcov(7,7)=1^2;
p_store=zeros(m,7);
p_store(1,:)=diag(pcov)';

% Main EKF Loop
for i=1:m-1

% Display When Every 500th Point Is Reached
if (i500==500), 
 disp(sprintf('      Filter has reached point %5i',i-1))
 i500=0;
end
i500=i500+1;    
    
% Update (we're updating at 0.1 seconds)
if kk == kup,

v1e=xe(i,1);v3e=xe(i,2);
alpe=atan(v3e/v1e);
w2e=xe(i,3);
thetae=xe(i,4);

dalpdv1=-v3e/(v1e^2+v3e^2);
dalpdv3=v1e/(v1e^2+v3e^2);

% Form H matrix
h=zeros(4,7);
h(1,1)=dalpdv1;
h(1,2)=dalpdv3;
h(2,1)=v1e*(v1e^2+v3e^2)^(-0.5);
h(2,2)=v3e*(v1e^2+v3e^2)^(-0.5);
h(3,3)=1;
h(4,4)=1;

% Estimated Output
ye=[atan(v3e/v1e);norm([v1e;v3e]);w2e;thetae];

% Gain
gain=pcov*h'*inv(h*pcov*h'+rcov);

% Covariance Update
pcov=(eye(7)-gain*h)*pcov;

% State Update
xe(i,:)=xe(i,:)+(gain*(ym(i,:)'-ye))';
kk=1;

end

kk=kk+1;

% Compute Variables for Partials
v1e=xe(i,1);v3e=xe(i,2);
alpe=atan(v3e/v1e);
w2e=xe(i,3);
thetae=xe(i,4);

dalpdv1=-v3e/(v1e^2+v3e^2);
dalpdv3=v1e/(v1e^2+v3e^2);
    
cd0e=xe(i,5);cl0e=xe(i,6);cm0e=xe(i,7);
dragc=cd0e+cda*alpe+cdde*de(i);
liftc=cl0e+cla*alpe+clde*de(i);
momc=cm0e+cma*alpe+cmde*de(i)+cmq*(w2e-w2_ss)*cbar/2/v_ss;

drage=0.5*dragc*rho*(v1e^2+v3e^2)*s;
lifte=0.5*liftc*rho*(v1e^2+v3e^2)*s;

% Partials
f=zeros(7);

ddragdv1=rho*dragc*v1e*s-rho/2*cda*v3e*s;
ddragdv3=rho*dragc*v3e*s+rho/2*cda*v1e*s;
dliftdv1=rho*liftc*v1e*s-rho/2*cla*v3e*s;
dliftdv3=rho*liftc*v3e*s+rho/2*cla*v1e*s;

f(1,1)=-1/mass*(thrust(i)-drage)*sin(alpe)*dalpdv1-1/mass*ddragdv1*cos(alpe)...
       +1/mass*lifte*cos(alpe)*dalpdv1+1/mass*dliftdv1*sin(alpe);

f(1,2)=-1/mass*(thrust(i)-drage)*sin(alpe)*dalpdv3-1/mass*ddragdv3*cos(alpe)...
       +1/mass*lifte*cos(alpe)*dalpdv3+1/mass*dliftdv3*sin(alpe)-w2e;
   
f(1,3)=-v3e;
f(1,4)=-g*cos(thetae);
f(1,5)=-1/2/mass*rho*(v1e^2+v3e^2)*s*cos(alpe);
f(1,6)=1/2/mass*rho*(v1e^2+v3e^2)*s*sin(alpe);

f(2,1)=1/mass*(thrust(i)-drage)*cos(alpe)*dalpdv1-1/mass*ddragdv1*sin(alpe)...
      +1/mass*lifte*sin(alpe)*dalpdv1-1/mass*dliftdv1*cos(alpe)+w2e;

f(2,2)=1/mass*(thrust(i)-drage)*cos(alpe)*dalpdv3-1/mass*ddragdv3*sin(alpe)...
      +1/mass*lifte*sin(alpe)*dalpdv3-1/mass*dliftdv3*cos(alpe);
 
f(2,3)=v1e;
f(2,4)=-g*sin(thetae);
f(2,5)=-1/2/mass*rho*(v1e^2+v3e^2)*s*sin(alpe);
f(2,6)=-1/2/mass*rho*(v1e^2+v3e^2)*s*cos(alpe);

f(3,1)=rho*s*cbar/in(2,2)*(momc*v1e-1/2*cma*v3e);
     
f(3,2)=rho*s*cbar/in(2,2)*(momc*v3e+1/2*cma*v1e);

f(3,3)=1/in(2,2)/4/v_ss*cmq*rho*s*cbar^2*(v1e^2+v3e^2);

f(3,7)=1/in(2,2)*rho*(v1e^2+v3e^2)*s*cbar;

f(4,3)=1;

% Covariance Propagation
q=zeros(7);
xx=[pcov(:,1);pcov(:,2);pcov(:,3);pcov(:,4);pcov(:,5);pcov(:,6);pcov(:,7)];
f1=dt*air_covfuns(xx,f,q);
f2=dt*air_covfuns(xx+0.5*f1,f,q);
f3=dt*air_covfuns(xx+0.25*f1+0.25*f2,f,q);
f4=dt*air_covfuns(xx-f2+2*f3,f,q);
f5=dt*air_covfuns(xx+7/27*f1+10/27*f2+1/27*f4,f,q);
f6=dt*air_covfuns(xx+28/625*f1-1/5*f2+546/625*f3+54/625*f4-378/625*f5,f,q);
xx=xx+(1/24*f1+5/48*f4+27/56*f5+125/336*f6);
pcov=[xx(1:7) xx(8:14) xx(15:21) xx(22:28) xx(29:35) xx(36:42) xx(43:49)];
p_store(i+1,:)=diag(pcov)';

% State Propagation
coef=[cd0e;cda;cdde;cl0e;cla;clde;cm0e;cma;cmde;cmq];
f1=dt*air_longfuns_ekf(xe(i,:),de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
f2=dt*air_longfuns_ekf(xe(i,:)+0.5*f1',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
f3=dt*air_longfuns_ekf(xe(i,:)+0.25*f1'+0.25*f2',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
f4=dt*air_longfuns_ekf(xe(i,:)-f2'+2*f3',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
f5=dt*air_longfuns_ekf(xe(i,:)+7/27*f1'+10/27*f2'+1/27*f4',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
f6=dt*air_longfuns_ekf(xe(i,:)+28/625*f1'-1/5*f2'+546/625*f3'+54/625*f4'-378/625*f5',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
xe(i+1,:)=xe(i,:)+(1/24*f1'+5/48*f4'+27/56*f5'+125/336*f6');

end

% 3-Sigma Outlier
sig3=p_store.^(0.5)*3;

k=[1:kup:m]';

% Plot Results
subplot(221)
plot(t(k),xe(k,5))
axis([0 30 -0.04 0.04])
set(gca,'fontsize',12)
set(gca,'xtick',[0 5 10 15 20 25 30]);
set(gca,'ytick',[-0.04 -0.02 0 0.02 0.04]);
ylabel('{\it C_{D_0}} Estimate')
xlabel('Time (Sec)')
grid

subplot(222)
plot(t(k),xe(k,6))
axis([0 30 0 0.3])
set(gca,'fontsize',12)
set(gca,'xtick',[0 5 10 15 20 25 30]);
set(gca,'ytick',[-0.02 0 0.1 0.2 0.3]);
ylabel('{\it C_{L_0}} Estimate')
xlabel('Time (Sec)')
grid

subplot(223)
plot(t(k),xe(k,7))
axis([0 30 -0.02 0.02])
set(gca,'fontsize',12)
set(gca,'xtick',[0 5 10 15 20 25 30]);
set(gca,'ytick',[-0.02 -0.01 0 0.01 0.02]);
ylabel('{\it C_{m_0}} Estimate')
xlabel('Time (Sec)')
grid

subplot(224)
plot(t(k),[sig3(k,4)*180/pi (xe(k,4)-x(k,4))*180/pi -sig3(k,4)*180/pi])
axis([0 30 -0.1 0.1])
set(gca,'fontsize',12)
set(gca,'xtick',[0 5 10 15 20 25 30]);
set(gca,'ytick',[-0.1 -0.05 0 0.05 0.1]);
ylabel('Pitch Error (Deg)')
xlabel('Time (Sec)')
grid