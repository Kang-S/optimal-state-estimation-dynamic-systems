function x=orbit_initial(y,tm,lat,sid);

% x=orbit_initial(y,tm,lat,sid);
%
% Determines the initial orbit conditions at tm(2) using
% a Herrick-Gibbs approach. The inputs are:
%     y = [range azimuth elevation] in km and rad (3x3)
%    tm = measurement times (assumed multiples of dt, and tm(1) = t0) (3x1)
%   lat = geocentric lattitude of radar in rad (1x1)
%   sid = sidereal time in sec (3x1)

% Written by John L. Crassidis 9/03

% Constants
mu=398600.64;
rearth=6378;
x=zeros(6,1);

% Compute Vectors
bigr1=rearth*[cos(lat)*cos(sid(1));cos(lat)*sin(sid(1));sin(lat)];
bigr2=rearth*[cos(lat)*cos(sid(2));cos(lat)*sin(sid(2));sin(lat)];
bigr3=rearth*[cos(lat)*cos(sid(3));cos(lat)*sin(sid(3));sin(lat)];

% Matrix Conversion
mc=[cos(lat) 0 sin(lat);0 1 0;-sin(lat) 0 cos(lat)];
m1=[cos(sid(1)) sin(sid(1)) 0;-sin(sid(1)) cos(sid(1)) 0; 0 0 1];
m2=[cos(sid(2)) sin(sid(2)) 0;-sin(sid(2)) cos(sid(2)) 0; 0 0 1];
m3=[cos(sid(3)) sin(sid(3)) 0;-sin(sid(3)) cos(sid(3)) 0; 0 0 1];

% Get rho_up, rho_east, rho_north
rho1=y(1,1)*[sin(y(1,3));cos(y(1,3))*sin(y(1,2));cos(y(1,3))*cos(y(1,2))];
rho2=y(2,1)*[sin(y(2,3));cos(y(2,3))*sin(y(2,2));cos(y(2,3))*cos(y(2,2))];
rho3=y(3,1)*[sin(y(3,3));cos(y(3,3))*sin(y(3,2));cos(y(3,3))*cos(y(3,2))];

% Compute Position Vector
r1=bigr1+(mc*m1)'*rho1;
r2=bigr2+(mc*m2)'*rho2;
r3=bigr3+(mc*m3)'*rho3;

% Modified Times (with k=1)
t12=(tm(2)-tm(1));
t13=(tm(3)-tm(1));
t23=(tm(3)-tm(2));
g1=t23/t12/t13;g3=t12/t23/t13;g2=g1-g3;
h1=mu*t23/12;h3=mu*t12/12;h2=h1-h3;

% Form Coefficents
d1=g1+h1/norm(r1)^3;
d2=g2+h2/norm(r2)^3;
d3=g3+h3/norm(r3)^3;

% Velocity
rd2=-d1*r1+d2*r2+d3*r3;

% Final Vector
x(1)=r2(1);x(2)=r2(2);x(3)=r2(3);
x(4)=rd2(1);x(5)=rd2(2);x(6)=rd2(3);