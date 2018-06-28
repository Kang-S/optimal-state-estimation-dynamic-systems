function f=air_fun(x,de,dr,da,thrust,coef,other,in,w1_ss,w2_ss,w3_ss,v_ss);

% Written by John L. Crassidis 9/03

% Function Routine for General Aircraft Equations

% Inertias
ixx=in(1,1);iyy=in(2,2);izz=in(3,3);ixz=in(1,3);

% Velocities and Euler Angles
v1=x(1);v2=x(2);v3=x(3);
w1=x(4);w2=x(5);w3=x(6);
phi=x(10);theta=x(11);psi=x(12);

% Get Other Constants
f=zeros(12,1);
m=other(5);
g=other(6);

% Speed, Angle of Attack and Sideslip
vm=norm([x(1);x(2);x(3)]);
alp=atan(x(3)/x(1));
bet=asin(x(2)/vm);

% Dynamic Pressure
q=0.5*other(1)*vm^2;

% General Coefficients
cd=coef(1)+coef(2)*alp+coef(3)*de;
cy=coef(4)+coef(5)*bet+coef(6)*dr+coef(7)*da;
cl=coef(8)+coef(9)*alp+coef(10)*de;

dd=2*vm^2;
cll=coef(11)+coef(12)*bet+coef(13)*dr+coef(14)*da+coef(23)*(w1-w1_ss)*other(4)/2/v_ss...
    +coef(24)*(w3-w3_ss)*other(4)/2/v_ss;
cm=coef(15)+coef(16)*alp+coef(17)*de+coef(22)*(w2-w2_ss)*other(3)/2/v_ss;
cn=coef(18)+coef(19)*bet+coef(20)*dr+coef(21)*da+coef(25)*(w1-w1_ss)*other(4)/2/v_ss...
    +coef(26)*w3*other(4)/2/v_ss;

% Drag, Force and Lift 
drag=cd*q*other(2);
yforce=cy*q*other(2);
lift=cl*q*other(2);

% Torques
la1=cll*q*other(2)*other(4);
la2=cm*q*other(2)*other(3);
la3=cn*q*other(2)*other(4);

% Forces
force1=-drag*cos(alp)+lift*sin(alp)+thrust*cos(alp);
force2=yforce;
force3=-drag*sin(alp)-lift*cos(alp)+thrust*sin(alp);

% Functions
f(1)=-g*sin(theta)+force1/m+v2*w3-v3*w2;
f(2)=g*cos(theta)*sin(phi)+force2/m+v3*w1-v1*w3;
f(3)=g*cos(theta)*cos(phi)+force3/m+v1*w2-v2*w1;

k1=ixz*(iyy-ixx);
k2=izz*(izz-iyy);
k3=ixz*(izz-iyy);
k4=ixx*(iyy-ixx);

f(4)=(ixz*la3+izz*la1-k1*w1*w2-ixz^2*w2*w3+ixz*izz*w1*w2-k2*w2*w3)/(ixx*izz-ixz^2);
f(5)=(la2-(ixx-izz)*w1*w2-ixz*(w1^2-w3^2))/iyy;
f(6)=(ixx*la3+ixz*la1+ixz^2*w1*w2-k3*w2*w3-k4*w1*w2-ixx*ixz*w2*w3)/(izz*ixx-ixz^2);

f(7)=cos(theta)*cos(psi)*v1...
     +(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi))*v2...
     +(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi))*v3;
f(8)=cos(theta)*sin(psi)*v1...
     +(sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi))*v2...
     +(cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi))*v3;
f(9)=-sin(theta)*v1+sin(phi)*cos(theta)*v2+cos(phi)*cos(theta)*v3;

f(10)=w1+sin(phi)*tan(theta)*w2+cos(phi)*tan(theta)*w3;
f(11)=cos(phi)*w2-sin(phi)*w3;
f(12)=sin(phi)*sec(theta)*w2+cos(phi)*sec(theta)*w3;