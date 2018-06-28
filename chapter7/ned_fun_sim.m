function f=ned_fun(x,acc_ned,wg);

% John L. Crassidis 3/04

% Pre-allocate
f=zeros(10,1);

% Semimajor Axis, Eccentricity and Earth Rate
a=6378137;
e=0.0818;
we=7.292155e-5;

% Get Variables
q=x(1:4)';
lat=x(5);
long=x(6);
height=x(7);
vn=x(8);
ve=x(9);
vd=x(10);

sin_lat=sin(lat);
cos_lat=cos(lat);

r_lat=a*(1-e^2)/((1-e^2*sin_lat^2)^(1.5));
r_long=a/((1-e^2*sin_lat^2)^(0.5));

den_lat=r_lat+height;
den_long=r_long+height;

% Gravity
g=9.780327*(1+0.0053024*sin_lat^2-0.0000058*sin(2*lat)^2)...
    -(3.0877e-6-0.0044e-6*sin_lat^2)*height+0.0072e-12*height^2;

% Quaternion
a_ned2body=attm(q);

w_body=wg(:)-a_ned2body*(we*[cos_lat;0;-sin_lat]+[ve/den_long;-vn/den_lat;-ve*sin_lat/den_long/cos_lat]);
om=[-crossm(w_body) w_body;-w_body' 0];
f(1:4)=0.5*om*q;

% Latitude, Longitude and Height Equations
f(5)=vn/den_lat;
f(6)=ve/den_long/cos_lat;
f(7)=-vd;

% Velocity Equations
f_acc=acc_ned(:);

f(8)=-(ve/den_long/cos_lat+2*we)*ve*sin_lat+vn*vd/den_lat+f_acc(1);
f(9)=(ve/den_long/cos_lat+2*we)*vn*sin_lat+ve*vd/den_long+2*we*vd*cos_lat+f_acc(2);
f(10)=-ve^2/den_long-vn^2/den_lat-2*we*ve*cos_lat+g+f_acc(3);