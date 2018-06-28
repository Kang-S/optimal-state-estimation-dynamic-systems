% Written by John L. Crassidis 9/03

% 747 Data, Used for Trim Valuse
cd00=0.0164;cda0=0.20;cdde0=0;
cl00=0.21;cla0=4.4;clde0=0.32;
cm00=0;cma0=-1.00;cmde0=-1.30;cmq0=-20.5;

% 747 Data in SI Units (Watch Velocity Below!!)
rho=.6536033;s=510.97;cbar=8.321;l=59.74;mass=2831897.6/9.81;
in=diag([24675882 44877565 67384138]);in(1,3)=1315143.1;in(3,1)=in(1,3);
g=9.81;

coef=[cd0;cda;cdde;cl0;cla;clde;cm0;cma;cmde;cmq];

other=[rho;s;cbar;l;mass;g];

% Initial Velocity
v_mag=205.13;

% Get Trim Values
qtrim=0.5*rho*v_mag^2;
dtrim=cla0*cmde0-cma0*clde0;
alptrim=((mass*g/qtrim/s-cl00)*cmde0+cm00*clde0)/dtrim;
detrim=(-cla0*cm00-cma0*(mass*g/qtrim/s-cl00))/dtrim;
dragtrim=(cd00+cda0*alptrim+cdde0*detrim)*qtrim*s;

% Initial Velocities [note: these must satisfy Eqs. (3.214a) and (3.215)]
v10=v_mag/sqrt(1+tan(alptrim)^2);
v30=v10*tan(alptrim);

% Other Initial Conditions
w20=0;
xx0=0;zz0=6096;

% Initial Pitch
theta0=alptrim;

% Steady-State Values
w2_ss=w20;
v_ss=v_mag;

% Sampling Interval
dt=0.2;
t=[0:dt:100]';
m=length(t);
x=zeros(m,6);
x(1,:)=[v10 v30 w20 xx0 zz0 theta0];

% Trim Conditions
de=detrim*ones(m,1);
thrust=dragtrim*ones(m,1);

de=detrim*ones(m,1);de(1:50)=detrim-ones(50,1)*1*pi/180;

% Main Loop
for i=1:m-1,
 f1=dt*air_longfun(x(i,:),de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
 f2=dt*air_longfun(x(i,:)+0.5*f1',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
 f3=dt*air_longfun(x(i,:)+0.5*f2',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
 f4=dt*air_longfun(x(i,:)+f3',de(i),thrust(i),coef,other,in(2,2),w2_ss,v_ss);
 x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
end

vel=x(:,1:2);
velmag=(vel(:,1).^2+vel(:,2).^2).^(0.5);
alp=atan(vel(:,2)./vel(:,1));

w=x(:,3);
pos=x(:,4:5);
theta=x(:,6);