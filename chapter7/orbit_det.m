function [xe,xecov]=orbit_det(x0,t0,tf,dt,y,tm,lat,sid,max,met,ymcov);

% [xe,xecov]=orbit_det(x0,t0,tf,dt,y,tm,rad_az0,rad_el,max,met,ymcov);
%
% Determines the initial orbit conditions at Epoch using
% a differential least-squares approach. The inputs are:
%    x0 = initial guess (6x1)
%    t0 = initial time (1x1)
%    tf = final time (1x1)
%    dt = integration interval in sec. (1x1)
%     y = [range azimuth elevation] in km and rad (mx3)
%    tm = measurement times (assumed multiples of dt, and tm(1) = t0) (mx3)
%   lat = geocentric lattitude of radar in rad (1x1)
%   sid = sidereal time in sec (mx1)
%   max = maximum number of iterations
%   met = 1 for discrete phi (fastest)
%       = 2 for Battin's method
%       = 3 for Escobol's method
% ymcov = measurement covariance (3x3)

% Written by John L. Crassidis 9/03

% Find Measurement Times
t=[t0:dt:tf]';
clear k
for i=1:length(tm),
 k(i)=find(t==tm(i));
end

% Range, Azimuth, and Elevation
ym=zeros(3*length(y),1);
ym(:)=[y(:,1) y(:,2) y(:,3)]';
w=inv(ymcov);
ww=kron(eye(length(k)),w);

% Initialize Matrices
ye1=zeros(length(y),1);
ye2=zeros(length(y),1);
ye3=zeros(length(y),1);
xe=x0(:);
phi=eye(6);
icount=0;xee=zeros(6,1);clear xev;xev(1,:)=xe';
phia=lat;theta=sid;
h=[];
rearth=6378;
mu=398600.64;

% Main loop
while norm(xe-xee)>.00000001

icount=icount+1;
if icount > max, break, end

% Propagator and State Transition Matrix
[bigphi,xstate]=orb_prop(xe,eye(6),t0,tf,dt,tm(2:length(tm)),met);

for jjj=1:length(tm);
 iii=find(t==tm(jjj));
 xw=xstate(iii,:);
if jjj==1,
 phi=eye(6);
else

if met==2 | met==3,
 phi=bigphi(6*(jjj-1)-5:6*(jjj-1),:);
else
 rw=norm(xw(1:3));
 asub=mu/rw^5*[3*xw(1)^2-rw^2 3*xw(1)*xw(2) 3*xw(1)*xw(3)
                  3*xw(1)*xw(2) 3*xw(2)^2-rw^2 3*xw(2)*xw(3)
                  3*xw(1)*xw(3 ) 3*xw(2)*xw(3) 3*xw(3)^2-rw^2]; 
 biga=[zeros(3) eye(3);asub zeros(3)];
 bigad=c2d(biga,zeros(6,1),tm(2)-tm(1));
 phi=bigad*phi;
end
end

% Estimated Quantities
rhoud=cos(phia)*cos(theta(jjj)).*(xw(1)-rearth*cos(phia)*cos(theta(jjj))) ...
    +cos(phia)*sin(theta(jjj)).*(xw(2)-rearth*cos(phia)*sin(theta(jjj))) ...
    +sin(phia)*(xw(3)-rearth*sin(phia));
rhoed=-sin(theta(jjj)).*(xw(1)-rearth*cos(phia)*cos(theta(jjj))) ...
     +cos(theta(jjj)).*(xw(2)-rearth*cos(phia)*sin(theta(jjj)));
rhond=-sin(phia)*cos(theta(jjj)).*(xw(1)-rearth*cos(phia)*cos(theta(jjj))) ...
     -sin(phia)*sin(theta(jjj)).*(xw(2)-rearth*cos(phia)*sin(theta(jjj))) ...
     +cos(phia)*(xw(3)-rearth*sin(phia));
rhod=norm([rhoud rhoed rhond]);

% Partials
drdx=(rhoud*cos(phia)*cos(theta(jjj))-rhoed*sin(theta(jjj)) ...
    -rhond*sin(phia)*cos(theta(jjj)))/rhod;

drdy=(rhoud*cos(phia)*sin(theta(jjj))+rhoed*cos(theta(jjj)) ...
    -rhond*sin(phia)*sin(theta(jjj)))/rhod;

drdz=(rhoud*sin(phia)+rhond*cos(phia))/rhod;


fac1=1/(1+(rhoed/rhond)^2);
dadx=fac1*(-sin(theta(jjj))+rhoed*sin(phia)*cos(theta(jjj))/rhond)/rhond;
dady=fac1*(cos(theta(jjj))+rhoed*sin(phia)*sin(theta(jjj))/rhond)/rhond;
dadz=fac1*(-rhoed*cos(phia)/rhond)/rhond;

fac2=1/sqrt(1-(rhoud/rhod)^2);
dedx=fac2*(cos(phia)*cos(theta(jjj))-rhoud*drdx/rhod)/rhod;
dedy=fac2*(cos(phia)*sin(theta(jjj))-rhoud*drdy/rhod)/rhod;
dedz=fac2*(sin(phia)-rhoud*drdz/rhod)/rhod;

% Sensitivity Matrix
hpre=[drdx drdy drdz 0 0 0
      dadx dady dadz 0 0 0
      dedx dedy dedz 0 0 0];


h=[h;hpre*phi];
ye1(jjj)=rhod;
ye2(jjj)=atan2(rhoed,rhond);
if (ye2(jjj)-y(jjj,2) > pi), ye2(jjj)=ye2(jjj)-2*pi; end
if (y(jjj,2)-ye2(jjj) > pi), ye2(jjj)=ye2(jjj)+2*pi; end
ye3(jjj)=asin(rhoud/rhod);

end

% Update with Least Squares
ye=zeros(3*length(y),1);
ye(:)=[ye1 ye2 ye3]';

xee=xe;
hh=h;
xe=xe+inv(h'*ww*h)*h'*ww*(ym-ye);
h=[];
xe'

xev(icount+1,:)=xe';

end

xe=xev;
xecov=inv(hh'*ww*hh);