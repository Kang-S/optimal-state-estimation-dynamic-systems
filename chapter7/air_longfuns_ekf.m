function f=air_longfuns_ekf(x,de,thrust,coef,other,in,w2_ss,v_ss);

% Written by John L. Crassidis 9/03

% Function Routine for Longitudinal Aircraft Equations in EKF

% Inertia
iyy=in;

% Velocities and Pitch Angle
v1=x(1);v3=x(2);w2=x(3);
theta=x(4);

% Get Other Constants
f=zeros(7,1);
m=other(5);
g=other(6);

% Speed and Angle of Attack
vm=norm([x(1);x(2)]);
alp=atan(x(2)/x(1));

% Dynamic Pressure
q=0.5*other(1)*vm^2;

% Drag, Lift and Moment Coefficients
cd=x(5)+coef(2)*alp+coef(3)*de;
cl=x(6)+coef(5)*alp+coef(6)*de;
cm=x(7)+coef(8)*alp+coef(9)*de+coef(10)*(w2-w2_ss)*other(3)/2/v_ss;

% Drag and Lift
drag=cd*q*other(2);
lift=cl*q*other(2);

% Torque
la2=cm*q*other(2)*other(3);

% Forces
force1=-drag*cos(alp)+lift*sin(alp)+thrust*cos(alp);
force3=-drag*sin(alp)-lift*cos(alp)+thrust*sin(alp);

% Functions
f(1)=-g*sin(theta)+force1/m-v3*w2;
f(2)=g*cos(theta)+force3/m+v1*w2;


f(3)=(la2)/iyy;

f(4)=w2;

f(5)=0;
f(6)=0;
f(7)=0;