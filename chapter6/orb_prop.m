function [phimf,xstate]=orb_prop(xi,phim0,t0,tf,dt,tm,met)

% [xe,xecov]=orb_prop(xi,phim0,t0,tf,dt,tm,met);
%
% Orbit and State Transition Matrix Propagation Routine
% The inputs are:
%    xi = initial state (6x1)
% phim0 = initial state matrix (6x6)
%    t0 = initial time (1x1)
%    tf = final time (1x1)
%    dt = integration interval in sec. (1x1)
%    tm = measurement times (assumed multiples of dt, and tm(1) = t0) (mx3)
%   met = 1 for discrete phi (fastest)
%       = 2 for Battin's method
%       = 3 for Escobol's method

% Written by John L. Crassidis 9/03

% Initialize Variables
t=t0;
time=[t0:dt:tf];
len=length(time);
lentm=length(tm);
xstate=zeros(len,6);xstate(1,:)=xi';
phim=phim0;
phimf=zeros(6*lentm,6);

mu=398600.64;
x0=xi(1);y0=xi(2);z0=xi(3);xd0=xi(4);yd0=xi(5);zd0=xi(6);
r0=sqrt(xi(1)^2+xi(2)^2+xi(3)^2);
rd0=(x0*xd0+y0*yd0+z0*zd0)/r0;
v0=sqrt(xi(4)^2+xi(5)^2+xi(6)^2);
a=inv(2/r0-v0^2/mu);


% Propagate Orbit
for i=1:len
 f1=dt*orbitfun(xstate(i,:),mu);
 f2=dt*orbitfun(xstate(i,:)+0.5*f1',mu);
 f3=dt*orbitfun(xstate(i,:)+0.5*f2',mu);
 f4=dt*orbitfun(xstate(i,:)+f3',mu);
 xstate(i+1,:)=xstate(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
end

% State Transition Propagation
for jjkk=1:lentm,
iii=find(time==tm(jjkk));
xf=xstate(iii,:)';

dadx0=2*a*a*x0/r0^3;dady0=2*a*a*y0/r0^3;dadz0=2*a*a*z0/r0^3;
dadxd0=2*a*a*xd0/mu;dadyd0=2*a*a*yd0/mu;dadzd0=2*a*a*zd0/mu;

dcdx0=v0^2*x0/(mu*r0);dcdy0=v0^2*y0/(mu*r0);dcdz0=v0^2*z0/(mu*r0);
dcdxd0=2*r0*xd0/mu;dcdyd0=2*r0*yd0/mu;dcdzd0=2*r0*zd0/mu;

se=r0*rd0/sqrt(mu*a);
dsdx0=xd0/sqrt(mu*a)-a*se*x0/(r0^3);
dsdy0=yd0/sqrt(mu*a)-a*se*y0/(r0^3);
dsdz0=zd0/sqrt(mu*a)-a*se*z0/(r0^3);
dsdxd0=x0/sqrt(mu*a)-a*se*xd0/mu;
dsdyd0=y0/sqrt(mu*a)-a*se*yd0/mu;
dsdzd0=z0/sqrt(mu*a)-a*se*zd0/mu;

u0=[x0/r0;y0/r0;z0/r0];

duxdx0=1/r0-(x0^2)/(r0^3);
duxdy0=-x0*y0/(r0^3);
duxdz0=-x0*z0/(r0^3);
duydx0=-x0*y0/(r0^3);
duydy0=1/r0-(y0^2)/(r0^3);
duydz0=-y0*z0/(r0^3);
duzdx0=-x0*z0/(r0^3);
duzdy0=-y0*z0/(r0^3);
duzdz0=1/r0-(z0^2)/(r0^3);

s0=[r0*xd0-rd0*x0;r0*yd0-rd0*y0;r0*zd0-rd0*z0];

dsxdx0=rd0*(x0^2/(r0^2)-1);
dsxdy0=(y0*xd0-x0*yd0)/r0+x0*y0*rd0/(r0^2);
dsxdz0=(z0*xd0-x0*zd0)/r0+x0*z0*rd0/(r0^2);
dsydx0=(x0*yd0-y0*xd0)/r0+x0*y0*rd0/(r0^2);
dsydy0=rd0*(y0^2/(r0^2)-1);
dsydz0=(z0*yd0-y0*zd0)/r0+z0*y0*rd0/(r0^2);
dszdx0=(x0*zd0-z0*xd0)/r0+x0*z0*rd0/(r0^2);
dszdy0=(y0*zd0-z0*yd0)/r0+y0*z0*rd0/(r0^2);
dszdz0=rd0*(z0^2/(r0^2)-1);

dsxdxd0=r0-x0^2/r0;dsxdyd0=-x0*y0/r0;dsxdzd0=-x0*z0/r0;
dsydxd0=-x0*y0/r0;dsydyd0=r0-y0^2/r0;dsydzd0=-y0*z0/r0;
dszdxd0=-x0*z0/r0;dszdyd0=-y0*z0/r0;dszdzd0=r0-z0^2/r0;

r=xf(1:3);rd=xf(4:6);x=r(1);y=r(2);z=r(3);rm=norm(r);
cosv_v0=r'*[x0;y0;z0]/(norm(r)*norm([x0;y0;z0]));
s=(x0*y-x*y0)/abs(x0*y-x*y0);
sinv_v0=s*sqrt(1-cosv_v0^2);

ce=1-r0/a;
p=a*(1-ce^2-se^2);
e_e0=asin(rm/sqrt(a*p)*sinv_v0-rm/p*(1-cosv_v0)*se);
m_m0=e_e0+se*(1-cos(e_e0))-ce*sin(e_e0);

a1=rm/a*cosv_v0+3/2*a*m_m0*(1/sqrt(a*p)*(p/r0-ce)*sinv_v0 ...
  +se*(1/p-1/r0)*(1-cosv_v0)-se/rm);
a2=1/sqrt(mu*p)*(rm/(2*a)*sinv_v0-3/2*sqrt(a*p)/r0*m_m0 ...
  +3/2*sqrt(a/p)*m_m0*(1-cosv_v0));
drda=a1*u0+a2*s0;


c1=-a+a*a*se/sqrt(a*p)*sinv_v0...
  +((rm*r0-a*a*se^2+2*a*rm*ce)/p-a*rm/r0)*(1-cosv_v0) ...
  +a*rm*(ce/p-1/r0)*sinv_v0^2 ...
  +a*a*rm*se/sqrt(a*p)*(2/r0-1/p*(1+ce))*(1-cosv_v0)*sinv_v0 ...
  +a*a*se*se*rm/p*(1/p-1/r0)*(1-cosv_v0)^2;
c2=a*rm/sqrt(mu*p)*(sinv_v0/r0+1/p*sqrt(a/p)*se*(1-cosv_v0)^2 ...
  -sinv_v0/p*(1-cosv_v0));
drdce=c1*u0+c2*s0;

s1=a*rm/sqrt(a*p)*sinv_v0*(1+(1-cosv_v0)*(1-r0/p*ce)) ...
  +se*a*rm/p*(1-cosv_v0)*cosv_v0 ...
  +se*a*r0/p*(1-cosv_v0)*(rm/p*(1-cosv_v0)-1);
s2=1/sqrt(mu*p)*(rm*r0*sqrt(a*p)/p/p*(1-cosv_v0)^2);
drdse=s1*u0+s2*s0;

xv=rm*cosv_v0;yv=rm*sinv_v0;

rdm=r'*rd/rm;
d=rm*rdm/sqrt(mu*a);
a1b=(3/2*m_m0/rm*sqrt(mu*p)*a*d*(1/r0-1/rm)-1/2*sqrt(mu*p)/a)*sinv_v0 ...
   -(3/2*m_m0/rm/rm*sqrt(mu*a)*a*d*d+3/2*m_m0/rm/r0*sqrt(mu*a)*a*d*se ...
   +3/2*m_m0*sqrt(mu/a)+1/2*sqrt(mu/a)*d)*(1-cosv_v0) ...
   +3/2*m_m0*sqrt(mu/a)*(1-a/rm)+3/2*m_m0/rm/r0*p*sqrt(mu*a)+1/2*sqrt(mu/a)*d;
a2b=(3/2*m_m0/a*(1-a/rm)*sqrt(a/p)+3/2*m_m0/rm/rm*sqrt(a/p)*a*d*d ...
   +3/2*m_m0/rm/r0*sqrt(a*p))*sinv_v0 ...
   +(3/2*m_m0/rm*a*d/p-3/2*m_m0/rm/r0*a*se-3/2*m_m0/rm/rm*a*d)*(1-cosv_v0) ...
   -3/2*m_m0/rm/r0*a*d;
drrdda=a1b*u0+a2b*s0;

beta1=sinv_v0/sqrt(a*p)-se/p*(1-cosv_v0);
beta2=a*(d*a/rm/sqrt(a*p)*sinv_v0-(1/r0+1/rm) ...
     +1/p*(r0/a-d*a*se/rm)*(1-cosv_v0)); 
beta3=(beta2+2*a*ce/p)*(1-cosv_v0)-p*a*a*beta1*beta1/r0;
beta4=-(a*ce/p+beta2)*sinv_v0-a*sqrt(a/p)*beta1*(1-cosv_v0);
c1b=(a*sqrt(mu*a)*beta1*cosv_v0+a*sqrt(mu/p)*ce*sinv_v0 ...
   +d*sqrt(mu*a)*beta3-sqrt(mu*p)*beta4);
c2b=(a*sqrt(a/p)*beta1*sinv_v0+d*a*a/p/sqrt(a*p)*ce*sinv_v0+beta3  ...
    +d*sqrt(a/p)*beta4);
drrddce=c1b*u0+c2b*s0;

gam1=sqrt(mu*a)*(1-r0/p*(1-cosv_v0));
gam2=a*(beta1-r0/rm/p*d*(1-cosv_v0));
gam3=(2*a*se/p+gam2+a*beta1)*(1-cosv_v0);
gam4=-(a*se/p+gam2)*sinv_v0+r0*sqrt(a*p)/p/p*(1-cosv_v0)^2;
s1b=gam1*cosv_v0+a*sqrt(mu/p)*se*sinv_v0+sqrt(mu*a)*d*gam3-sqrt(mu*p)*gam4;
s2b=gam1/sqrt(mu*p)*sinv_v0+a*a*se*d/p/sqrt(a*p)*sinv_v0+gam3+d*sqrt(a/p)*gam4;
drrddse=s1b*u0+s2b*s0;

ub=sqrt(mu*a)*d*cosv_v0-sqrt(mu*p)*sinv_v0;
drrddux0=[ub;0;0];drrdduy0=[0;ub;0];drrdduz0=[0;0;ub];

sb=cosv_v0+sqrt(a/p)*d*sinv_v0;
drrddsx0=[sb;0;0];drrddsy0=[0;sb;0];drrddsz0=[0;0;sb];

drmda=1/rm*(x*drda(1)+y*drda(2)+z*drda(3));
drdda=1/rm*(drrdda-rd*drmda);

drmdce=1/rm*(x*drdce(1)+y*drdce(2)+z*drdce(3));
drddce=1/rm*(drrddce-rd*drmdce);

drmdse=1/rm*(x*drdse(1)+y*drdse(2)+z*drdse(3));
drddse=1/rm*(drrddse-rd*drmdse);

drmdux0=1/rm*x*xv;
drddux0=1/rm*(drrddux0-rd*drmdux0);

drmduy0=1/rm*y*xv;
drdduy0=1/rm*(drrdduy0-rd*drmduy0);

drmduz0=1/rm*z*xv;
drdduz0=1/rm*(drrdduz0-rd*drmduz0);

drmdsx0=1/rm*x*yv/sqrt(mu*p);
drddsx0=1/rm*(drrddsx0-rd*drmdsx0);

drmdsy0=1/rm*y*yv/sqrt(mu*p);
drddsy0=1/rm*(drrddsy0-rd*drmdsy0);

drmdsz0=1/rm*z*yv/sqrt(mu*p);
drddsz0=1/rm*(drrddsz0-rd*drmdsz0);


xb=[drda(1) drdce(1) drdse(1) xv 0 0 yv/sqrt(mu*p) 0 0]';
yb=[drda(2) drdce(2) drdse(2) 0 xv 0 0 yv/sqrt(mu*p) 0]';
zb=[drda(3) drdce(3) drdse(3) 0 0 xv 0 0 yv/sqrt(mu*p)]';
xbd=[drdda(1) drddce(1) drddse(1) drddux0(1) drdduy0(1) drdduz0(1) ...
     drddsx0(1) drddsy0(1) drddsz0(1)]';
ybd=[drdda(2) drddce(2) drddse(2) drddux0(2) drdduy0(2) drdduz0(2) ...
     drddsx0(2) drddsy0(2) drddsz0(2)]';
zbd=[drdda(3) drddce(3) drddse(3) drddux0(3) drdduy0(3) drdduz0(3) ...
     drddsx0(3) drddsy0(3) drddsz0(3)]';

qx=[dadx0 dcdx0 dsdx0 duxdx0 duydx0 duzdx0 dsxdx0 dsydx0 dszdx0]';
qy=[dady0 dcdy0 dsdy0 duxdy0 duydy0 duzdy0 dsxdy0 dsydy0 dszdy0]';
qz=[dadz0 dcdz0 dsdz0 duxdz0 duydz0 duzdz0 dsxdz0 dsydz0 dszdz0]';
qxd=[dadxd0 dcdxd0 dsdxd0 0 0 0 dsxdxd0 dsydxd0 dszdxd0]';
qyd=[dadyd0 dcdyd0 dsdyd0 0 0 0 dsxdyd0 dsydyd0 dszdyd0]';
qzd=[dadzd0 dcdzd0 dsdzd0 0 0 0 dsxdzd0 dsydzd0 dszdzd0]';


bigm=[xb'*qx xb'*qy xb'*qz 
      yb'*qx yb'*qy yb'*qz 
      zb'*qx zb'*qy zb'*qz];
bign=[xb'*qxd xb'*qyd xb'*qzd
      yb'*qxd yb'*qyd yb'*qzd
      zb'*qxd zb'*qyd zb'*qzd];

bigms=[xbd'*qx xbd'*qy xbd'*qz 
       ybd'*qx ybd'*qy ybd'*qz 
       zbd'*qx zbd'*qy zbd'*qz];


bigns=[xbd'*qxd xbd'*qyd xbd'*qzd 
       ybd'*qxd ybd'*qyd ybd'*qzd 
       zbd'*qxd zbd'*qyd zbd'*qzd];

phif=[bigm bign;bigms bigns];


if met==3,
 phimf(6*jjkk-5:6*jjkk,1:6)=phif;
end

h=cross(r(:),rd(:));h=sqrt(h'*h);
p=h*h/mu;sig0=[x0 y0 z0]*[xd0;yd0;zd0]/sqrt(mu);

theta=acos(cosv_v0);

ccc=cross([xd0;yd0;zd0],cross(r(:),rd(:)))-mu*[x0;y0;z0]/r0;
ecc=norm(ccc)/mu;
ie=ccc/(mu*ecc);
theta=acos(r(:)'*ie/norm(r(:))/norm(ie))-acos([x0;y0;z0]'*ie/r0/norm(ie));

f=1-rm/p*(1-cos(theta));
g=rm*r0/sqrt(mu*p)*sin(theta);
ft=sqrt(mu)/r0/p*(sig0*(1-cos(theta))-sqrt(p)*sin(theta));
gt=1-r0/p*(1-cos(theta));

sig=r(:)'*rd(:)/sqrt(mu);
alp=1/a;

chi=alp*sqrt(mu)*(time(iii)-t0)+sig-sig0;

u2=(1-cos(sqrt(alp)*chi))/alp;
u3=(sqrt(alp)*chi-sin(sqrt(alp)*chi))/alp/sqrt(alp);
u4=chi^2/2/alp-u2/alp;
u5=chi^3/6/alp-u3/alp;

c=(3*u5-chi*u4-sqrt(mu)*(time(iii)-t0)*u2)/sqrt(mu);

r=r(:);rd=rd(:);
rv0=[x0;y0;z0];rdv0=[xd0;yd0;zd0];
bigr=r0/mu*(1-f)*((r-rv0)*rdv0'-(rd-rdv0)*rv0')+c/mu*rd*rdv0'+g*eye(3);
bigv=r0/mu*(rd-rdv0)*(rd-rdv0)'+1/rm^3*(r0*(1-f)*r*rv0'-c*r*rdv0')+gt*eye(3);
bigrb=rm/mu*(rd-rdv0)*(rd-rdv0)'+1/r0^3*(r0*(1-f)*r*rv0'+c*rd*rv0')+f*eye(3);
bigvb=-1/r0^2*(rd-rdv0)*rv0'-1/rm^2*r*(rd-rdv0)'...
      +ft*(eye(3)-1/rm^2*r*r'+1/mu/rm*(r*rd'-rd*r')*r*(rd-rdv0)') ...
      -mu*c/(rm^3*r0^3)*r*rv0';
phif=[bigrb bigr;bigvb bigv];

if met==2,
 phimf(6*jjkk-5:6*jjkk,1:6)=phif;
end

end