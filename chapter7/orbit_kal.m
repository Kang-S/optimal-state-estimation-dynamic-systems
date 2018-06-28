function [xe,xecov]=orbit_kal(x0,t0,tf,dt,y,tm,lat,sid,max,ymcov,p0,q);

% [xe,xecov]=orbit_kal(x0,t0,tf,dt,y,tm,lat,dis,max,ymcov,p0,q);
%
% Determines the initial orbit conditions at Epoch using
% an iterated Kalman filter approach. The inputs are:
%    x0 = initial guess (6x1)
%    t0 = initial time (1x1)
%    tf = final time (1x1)
%    dt = integration interval in sec. (1x1)
%     y = [range azimuth elevation] in km and rad (mx3)
%    tm = measurement times (assumed multiples of dt, and tm(1) = t0) (mx3)
%   lat = geocentric lattitude of radar in rad (1x1)
%   sid = sidereal time in sec (mx1)
%   max = maximum number of iterations
% ymcov = measurement covariance (3x3)
%    p0 = initial error covariance (6x6)
%     q = discrete process noise covariance (6x6)

% Written by John L. Crassidis 9/03

% Find Measurement Times
t=[t0:dt:tf]';
clear k
for i=1:length(tm),
k(i)=find(t==tm(i));
end
k=k(:);

% Measurements
ym=[y(:,1) y(:,2) y(:,3)];

% Initialize
m=length(t);
ye1=zeros(length(y),1);
ye2=zeros(length(y),1);
ye3=zeros(length(y),1);
mu=398600.64;
pcov=p0;
kkk=1;
xe=zeros(length(t),6);
xe(1,:)=x0(:)';
pf=zeros(length(t),6);
pb=zeros(length(t),6);
pf(1,:)=diag(pcov)';
clear xiter piter
rearth=6378;
phia=lat;theta=sid;

% Main Loop
for jjjj=1:max,

for i=1:m-1;

if (i==k(kkk)),

xw=xe(i,:);


% Estimated Quantities
rhoud=cos(phia)*cos(theta(kkk)).*(xw(1)-rearth*cos(phia)*cos(theta(kkk))) ...
    +cos(phia)*sin(theta(kkk)).*(xw(2)-rearth*cos(phia)*sin(theta(kkk))) ...
    +sin(phia)*(xw(3)-rearth*sin(phia));
rhoed=-sin(theta(kkk)).*(xw(1)-rearth*cos(phia)*cos(theta(kkk))) ...
     +cos(theta(kkk)).*(xw(2)-rearth*cos(phia)*sin(theta(kkk)));
rhond=-sin(phia)*cos(theta(kkk)).*(xw(1)-rearth*cos(phia)*cos(theta(kkk))) ...
     -sin(phia)*sin(theta(kkk)).*(xw(2)-rearth*cos(phia)*sin(theta(kkk))) ...
     +cos(phia)*(xw(3)-rearth*sin(phia));
rhod=norm([rhoud rhoed rhond]);

% Partials
drdx=(rhoud*cos(phia)*cos(theta(kkk))-rhoed*sin(theta(kkk)) ...
    -rhond*sin(phia)*cos(theta(kkk)))/rhod;

drdy=(rhoud*cos(phia)*sin(theta(kkk))+rhoed*cos(theta(kkk)) ...
    -rhond*sin(phia)*sin(theta(kkk)))/rhod;

drdz=(rhoud*sin(phia)+rhond*cos(phia))/rhod;
 

fac1=1/(1+(rhoed/rhond)^2);
dadx=fac1*(-sin(theta(kkk))+rhoed*sin(phia)*cos(theta(kkk))/rhond)/rhond;
dady=fac1*(cos(theta(kkk))+rhoed*sin(phia)*sin(theta(kkk))/rhond)/rhond;
dadz=fac1*(-rhoed*cos(phia)/rhond)/rhond;

fac2=1/sqrt(1-(rhoud/rhod)^2);
dedx=fac2*(cos(phia)*cos(theta(kkk))-rhoud*drdx/rhod)/rhod;
dedy=fac2*(cos(phia)*sin(theta(kkk))-rhoud*drdy/rhod)/rhod;
dedz=fac2*(sin(phia)-rhoud*drdz/rhod)/rhod;

% Sensitivity
hpre=[drdx drdy drdz 0 0 0
      dadx dady dadz 0 0 0
      dedx dedy dedz 0 0 0];

% Gain
gain=pcov*hpre'*inv(hpre*pcov*hpre'+ymcov);

% Estimates
ye1=rhod;
ye2=atan2(rhoed,rhond);
ye3=asin(rhoud/rhod);
if (ye2-y(kkk,2) > pi), ye2=ye2-2*pi; end
if (y(kkk,2)-ye2 > pi), ye2=ye2+2*pi; end

ye=[ye1 ye2 ye3]';

% Unpdate
 xe(i,:)=xe(i,:)+(gain*(ym(kkk,:)'-ye))';
 kkk=kkk+1;
 pcov=(eye(6)-gain*hpre)*pcov;

end

% Propagate
re=norm(xe(i,1:3));
mu=398600.64;
asub=mu/re^5*[3*xe(i,1)^2-re^2 3*xe(i,1)*xe(i,2) 3*xe(i,1)*xe(i,3)
                  3*xe(i,1)*xe(i,2) 3*xe(i,2)^2-re^2 3*xe(i,2)*xe(i,3)
                  3*xe(i,1)*xe(i,3) 3*xe(i,2)*xe(i,3) 3*xe(i,3)^2-re^2]; 
biga=[zeros(3) eye(3);asub zeros(3)];
ad=c2d(biga,zeros(6,1),dt);
pcov=ad*pcov*ad'+q;

% Integrate
f1=dt*orbitfun(xe(i,:),mu);
f2=dt*orbitfun(xe(i,:)+0.5*f1',mu);
f3=dt*orbitfun(xe(i,:)+0.5*f2',mu);
f4=dt*orbitfun(xe(i,:)+f3',mu);
xe(i+1,:)=xe(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');

pf(i+1,:)=diag(pcov)';

end

xf=xe;pb(m,:)=pf(m,:);pcov=p0;
kkk=kkk-1;

% Backwards
for i=m:-1:2,

if (i==k(kkk)),
 
% Estimate quantities
xw=xe(i,:);

rhoud=cos(phia)*cos(theta(kkk)).*(xw(1)-rearth*cos(phia)*cos(theta(kkk))) ...
    +cos(phia)*sin(theta(kkk)).*(xw(2)-rearth*cos(phia)*sin(theta(kkk))) ...
    +sin(phia)*(xw(3)-rearth*sin(phia));
rhoed=-sin(theta(kkk)).*(xw(1)-rearth*cos(phia)*cos(theta(kkk))) ...
     +cos(theta(kkk)).*(xw(2)-rearth*cos(phia)*sin(theta(kkk)));
rhond=-sin(phia)*cos(theta(kkk)).*(xw(1)-rearth*cos(phia)*cos(theta(kkk))) ...
     -sin(phia)*sin(theta(kkk)).*(xw(2)-rearth*cos(phia)*sin(theta(kkk))) ...
     +cos(phia)*(xw(3)-rearth*sin(phia));
rhod=norm([rhoud rhoed rhond]);

% Partial Derivatives
drdx=(rhoud*cos(phia)*cos(theta(kkk))-rhoed*sin(theta(kkk)) ...
    -rhond*sin(phia)*cos(theta(kkk)))/rhod;
drdy=(rhoud*cos(phia)*sin(theta(kkk))+rhoed*cos(theta(kkk)) ...
    -rhond*sin(phia)*sin(theta(kkk)))/rhod;

drdz=(rhoud*sin(phia)+rhond*cos(phia))/rhod;
 

fac1=1/(1+(rhoed/rhond)^2);
dadx=fac1*(-sin(theta(kkk))+rhoed*sin(phia)*cos(theta(kkk))/rhond)/rhond;
dady=fac1*(cos(theta(kkk))+rhoed*sin(phia)*sin(theta(kkk))/rhond)/rhond;
dadz=fac1*(-rhoed*cos(phia)/rhond)/rhond;

fac2=1/sqrt(1-(rhoud/rhod)^2);
dedx=fac2*(cos(phia)*cos(theta(kkk))-rhoud*drdx/rhod)/rhod;
dedy=fac2*(cos(phia)*sin(theta(kkk))-rhoud*drdy/rhod)/rhod;
dedz=fac2*(sin(phia)-rhoud*drdz/rhod)/rhod;

% Sensitivity
hpre=[drdx drdy drdz 0 0 0
      dadx dady dadz 0 0 0
      dedx dedy dedz 0 0 0];

% Estimates
ye1=rhod;
ye2=atan2(rhoed,rhond);
ye3=asin(rhoud/rhod);
if (ye2-y(kkk,2) > pi), ye2=ye2-2*pi; end
if (y(kkk,2)-ye2 > pi), ye2=ye2+2*pi; end

% Gain
gain=pcov*hpre'*inv(hpre*pcov*hpre'+ymcov);

% Update
ye=[ye1 ye2 ye3]';
xe(i,:)=xe(i,:)+(gain*(ym(kkk,:)'-ye))';
kkk=kkk-1;
pcov=(eye(6)-gain*hpre)*pcov;

end

% Propagate
re=norm(xe(i,1:3));
mu=398600.64;
asub=mu/re^(5)*[3*xe(i,1)^2-re^2 3*xe(i,1)*xe(i,2) 3*xe(i,1)*xe(i,3)
                  3*xe(i,1)*xe(i,2) 3*xe(i,2)^2-re^2 3*xe(i,2)*xe(i,3)
                  3*xe(i,1)*xe(i,3) 3*xe(i,2)*xe(i,3) 3*xe(i,3)^2-re^2]; 
biga=[zeros(3) eye(3);asub zeros(3)];

ad=c2d(-biga,zeros(6,1),dt);
pcov=ad*pcov*ad'+q;

pb(i-1,:)=diag(pcov)';

% Integrate
f1=-dt*orbitfun(xe(i,:),mu);
f2=-dt*orbitfun(xe(i,:)+0.5*f1',mu);
f3=-dt*orbitfun(xe(i,:)+0.5*f2',mu);
f4=-dt*orbitfun(xe(i,:)+f3',mu);
xe(i-1,:)=xe(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');

end

plast=pcov;
xb=xe;pcov=p0;

xiter(jjjj,:)=xe(1,:);
piter(jjjj,:)=pb(1,:);

xe(1,:)

end

xe=xiter;
xecov=plast;




